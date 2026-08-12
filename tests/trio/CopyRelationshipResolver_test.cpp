#include "src/trio/impl/CopyRelationshipResolver.h"

#include "src/trio/model/AlignmentSiteGraph.h"
#include "gtest/gtest.h"

#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

void addResolverPair(const std::string &a, const std::string &b,
                     const std::string &sourceA, const std::string &sourceB,
                     const std::string &textA, const std::string &textB,
                     const std::string &run, AlignmentEvidenceSet &evidence) {
    std::ostringstream maf;
    maf << "a score=10\n"
        << "s " << sourceA << " 0 " << textA.size() << " + " << textA.size()
        << ' ' << textA << "\n"
        << "s " << sourceB << " 0 " << textB.size() << " + " << textB.size()
        << ' ' << textB << "\n";
    std::istringstream stream(maf.str());
    PairwiseMafInput input;
    input.leftTaxon = a; input.rightTaxon = b; input.mafPath = run; input.runId = run;
    PairwiseEvidenceNormalizer::appendBlocks(
        input, MafEvidenceReader::read(stream, run, a, b), evidence);
}

CopyConstraint interval(const std::string &taxon, const std::string &sequence,
                        const std::string &group, ConstraintStrength strength) {
    CopyConstraint value;
    value.recordType = CopyRecordType::MEMBER;
    value.familyId = "F";
    value.ancestralCopyId = group;
    value.taxonId = taxon;
    value.memberType = CopyMemberType::INTERVAL;
    value.hasInterval = true;
    value.interval = GenomicInterval(sequence, 0, 100, '+');
    value.relation = CopyRelation::ORTHOLOG;
    value.confidence = 1.0;
    value.strength = strength;
    value.source = "test";
    value.sourceLocation = SourceLocation("copy.tsv", 1);
    return value;
}

CopyConstraint anchor(const std::string &taxon, const std::string &memberId) {
    CopyConstraint value;
    value.recordType = CopyRecordType::MEMBER;
    value.familyId = "F";
    value.ancestralCopyId = "A";
    value.taxonId = taxon;
    value.memberType = CopyMemberType::ANCHOR;
    value.memberId = memberId;
    value.relation = CopyRelation::ORTHOLOG;
    value.confidence = 1.0;
    value.strength = ConstraintStrength::HARD;
    value.source = "test";
    value.sourceLocation = SourceLocation("copy.tsv", 17);
    return value;
}

TEST(CopyRelationshipResolver, AppliesHardIntervalsBeforeGraphConstruction) {
    AlignmentEvidenceSet evidence;
    addResolverPair("I1", "I2", "a", "b", "AC", "AC", "12", evidence);
    addResolverPair("I1", "O1", "a", "o", "AC", "AC", "1o", evidence);
    addResolverPair("I2", "O1", "b", "o", "AC", "AC", "2o", evidence);
    CopyConstraintSet constraints;
    constraints.records = {interval("I1", "a", "A", ConstraintStrength::HARD),
                           interval("I2", "b", "A", ConstraintStrength::HARD),
                           interval("O1", "o", "A", ConstraintStrength::HARD)};
    const CopyResolutionResult resolved = CopyRelationshipResolver::resolve(
        evidence, constraints, CopyResolverOptions());
    EXPECT_EQ(evidence.residues.size(), resolved.assignments.size());
    EXPECT_EQ(evidence.homologies.size(), resolved.selectedEvidenceIds.size());
    for (const auto &assignment : resolved.assignments) {
        EXPECT_EQ("F:A", assignment.second.copyGroup);
        EXPECT_TRUE(assignment.second.hard);
    }

    AlignmentGraphBuildOptions build;
    build.ingroup1 = "I1"; build.ingroup2 = "I2"; build.primaryOutgroup = "O1";
    build.copyAssignments = resolved.assignments;
    build.allowedEvidenceIds = resolved.selectedEvidenceIds;
    build.restrictToAllowedEvidence = true;
    build.requireResolvedCopyAssignments = true;
    EXPECT_EQ(2u, AlignmentSiteGraphBuilder::build(evidence, build).sites().size());
}

TEST(CopyRelationshipResolver, RejectsOverlappingIncompatibleHardIntervals) {
    AlignmentEvidenceSet evidence;
    addResolverPair("I1", "I2", "a", "b", "A", "A", "12", evidence);
    CopyConstraintSet constraints;
    constraints.records = {interval("I1", "a", "A", ConstraintStrength::HARD),
                           interval("I1", "a", "B", ConstraintStrength::HARD)};
    EXPECT_THROW(CopyRelationshipResolver::resolve(
                     evidence, constraints, CopyResolverOptions()),
                 CopyResolutionError);
}

TEST(CopyRelationshipResolver, RejectsAnchorMembersUntilTheyAreExpandedToIntervals) {
    AlignmentEvidenceSet evidence;
    addResolverPair("I1", "I2", "a", "b", "A", "A", "12", evidence);
    CopyConstraintSet constraints;
    constraints.records.push_back(anchor("I1", "gene1"));
    try {
        CopyRelationshipResolver::resolve(evidence, constraints,
                                          CopyResolverOptions());
        FAIL() << "expected unexpanded anchor constraint to be rejected";
    } catch (const CopyResolutionError &error) {
        const std::string message(error.what());
        EXPECT_NE(std::string::npos, message.find("anchor MEMBER constraint 'gene1'"));
        EXPECT_NE(std::string::npos, message.find("copy.tsv:17"));
        EXPECT_NE(std::string::npos,
                  message.find("member_type=interval"));
    }
}

TEST(CopyRelationshipResolver, CopyCapacityAloneDoesNotSelectAParalog) {
    AlignmentEvidenceSet evidence;
    addResolverPair("I1", "I2", "a", "b", "A", "A", "12", evidence);
    CopyConstraintSet one;
    CopyCapacity capacity;
    capacity.key.taxonId = "I1"; capacity.count = 1;
    one.defaultCapacities[capacity.key] = capacity;
    CopyConstraintSet many = one;
    many.defaultCapacities.begin()->second.count = 8;
    const CopyResolutionResult a = CopyRelationshipResolver::resolve(
        evidence, one, CopyResolverOptions());
    const CopyResolutionResult b = CopyRelationshipResolver::resolve(
        evidence, many, CopyResolverOptions());
    ASSERT_EQ(a.assignments.size(), b.assignments.size());
    auto ai = a.assignments.begin();
    auto bi = b.assignments.begin();
    for (; ai != a.assignments.end(); ++ai, ++bi) {
        EXPECT_EQ(ai->second.copyGroup, bi->second.copyGroup);
        EXPECT_EQ("unambiguous_copy-compatible_pairwise_component",
                  ai->second.provenance);
    }
}

TEST(CopyRelationshipResolver, StrictModeRequiresExplicitHardMembership) {
    AlignmentEvidenceSet evidence;
    addResolverPair("I1", "I2", "a", "b", "A", "A", "12", evidence);
    CopyResolverOptions options;
    options.strictExplicitGroups = true;
    EXPECT_THROW(CopyRelationshipResolver::resolve(
                     evidence, CopyConstraintSet(), options),
                 CopyResolutionError);
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
