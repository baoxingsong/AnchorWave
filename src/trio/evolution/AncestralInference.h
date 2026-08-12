#pragma once

#include "src/trio/Phylogeny.h"
#include "src/trio/model/AlignmentSiteGraph.h"

#include <cstddef>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

struct EvolutionInferenceConfig {
    std::string targetNode;
    NucleotideParsimonyModel nucleotideModel;
    BinaryParsimonyModel presenceModel;
    bool requireCopyResolved = true;
    bool includeAmbiguousPresence = false;
};

struct AncestralSiteCall {
    std::string siteId;
    std::string copyGroup;
    NucleotideInferenceResult nucleotide;
    BinaryInferenceResult presence;
    bool emitted = false;
    char emittedBase = 'N';
    std::string decision;
};

struct AncestralBlock {
    std::string blockId;
    std::string copyGroup;
    std::vector<std::string> siteIds;
    std::string sequence;
};

struct AncestralOverlay {
    std::string immutableCoreHash;
    std::string targetNode;
    std::map<std::string, AncestralSiteCall> calls;
    std::vector<AncestralBlock> blocks;

    void validateAgainst(const AlignmentSiteGraph &graph) const;
};

struct ChildProjectionRow {
    std::string occurrencePath;
    std::string taxon;
    std::string source;
    int64_t start0 = 0;
    int64_t size = 0;
    int64_t sourceSize = 0;
    std::string text;
    bool emptyComponent = false;
};

struct AncestorChildAlignmentBlock {
    std::string projectionId;
    std::string ancestorBlockId;
    int64_t ancestorStart0 = 0;
    std::string ancestorText;
    std::vector<std::string> siteIds;
    std::vector<ChildProjectionRow> children;
};

class AncestralInferenceError : public std::runtime_error {
public:
    explicit AncestralInferenceError(const std::string &message)
        : std::runtime_error(message) {}
};

class AncestralInference {
public:
    static AncestralOverlay infer(const AlignmentSiteGraph &graph,
                                  const Phylogeny &tree,
                                  const EvolutionInferenceConfig &config);

    static std::vector<AncestorChildAlignmentBlock> projectToIngroupChildren(
        const AlignmentSiteGraph &graph, const AncestralOverlay &overlay);
};

class AncestralWriters {
public:
    static void writeCalls(const AncestralOverlay &overlay, std::ostream &output);
    static void writeFasta(const AncestralOverlay &overlay, std::ostream &output);
    static void writeChildMaf(
        const AncestralOverlay &overlay,
        const std::vector<AncestorChildAlignmentBlock> &projections,
        std::ostream &output);
    static void writeChildMap(
        const std::vector<AncestorChildAlignmentBlock> &projections,
        std::ostream &output);
};

}  // namespace trio
}  // namespace anchorwave

