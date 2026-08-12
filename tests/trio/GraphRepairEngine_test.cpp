#include "src/trio/impl/GraphRepairEngine.h"

#include "src/trio/impl/CopyRelationshipResolver.h"
#include "gtest/gtest.h"

#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

void addRepairPair(const std::string &a, const std::string &b,
                   const std::string &left, const std::string &right,
                   const std::string &run, AlignmentEvidenceSet &evidence) {
    std::size_t leftBases = 0, rightBases = 0;
    for (char c : left) leftBases += c != '-';
    for (char c : right) rightBases += c != '-';
    std::ostringstream maf;
    maf << "a score=1\n"
        << "s " << a << "chr 0 " << leftBases << " + 100 " << left << "\n"
        << "s " << b << "chr 0 " << rightBases << " + 100 " << right << "\n";
    std::istringstream stream(maf.str());
    PairwiseMafInput input;
    input.leftTaxon = a; input.rightTaxon = b; input.runId = run; input.mafPath = run;
    PairwiseEvidenceNormalizer::appendBlocks(
        input, MafEvidenceReader::read(stream, run, a, b), evidence);
}

TEST(GraphRepairEngine, ChangesOnlyBoundedDiscordantCore) {
    AlignmentEvidenceSet evidence;
    // Shared A/T flanks; I1-O1 shifts the middle G to a gap while I2-O1 aligns it.
    addRepairPair("I1", "I2", "ACGT", "ACGT", "12", evidence);
    addRepairPair("I1", "O1", "AC-GT", "ACG-T", "1o", evidence);
    addRepairPair("I2", "O1", "ACGT", "ACGT", "2o", evidence);
    const CopyResolutionResult copies = CopyRelationshipResolver::resolve(
        evidence, CopyConstraintSet(), CopyResolverOptions());
    AlignmentGraphBuildOptions build;
    build.ingroup1 = "I1"; build.ingroup2 = "I2"; build.primaryOutgroup = "O1";
    build.copyAssignments = copies.assignments;
    build.allowedEvidenceIds = copies.selectedEvidenceIds;
    build.restrictToAllowedEvidence = true;
    build.requireResolvedCopyAssignments = true;
    const AlignmentSiteGraph before = AlignmentSiteGraphBuilder::build(evidence, build);
    ASSERT_FALSE(before.conflicts().empty());
    GraphRepairOptions options;
    options.poa.minimumScoreDelta = -1000.0;
    const GraphRepairResult repaired =
        GraphRepairEngine::repairEligibleConflicts(before, options);
    ASSERT_FALSE(repaired.audit.empty());
    bool applied = false;
    for (const GraphRepairAudit &audit : repaired.audit) {
        if (audit.disposition == "localized_poa_applied") {
            applied = true;
            EXPECT_EQ(audit.outsideHashBefore, audit.outsideHashAfter);
            EXPECT_FALSE(audit.immutableLeftSite.empty());
            EXPECT_FALSE(audit.immutableRightSite.empty());
        }
    }
    EXPECT_TRUE(applied);
    repaired.graph.validateInvariants();
    for (const PathSegment &path : repaired.graph.pathSegments()) {
        EXPECT_EQ(static_cast<std::size_t>(path.end0 - path.start0),
                  path.sequenceText.size());
    }
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
