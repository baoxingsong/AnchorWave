#pragma once

#include "src/trio/Phylogeny.h"
#include "src/trio/model/AlignmentSiteGraph.h"

#include <cstdint>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

/**
 * A structure-only block definition.
 *
 * siteIds establishes the canonical orientation of the block.  These blocks
 * are deliberately independent of AncestralBlock: no ancestral base call or
 * inferred sequence is required to describe a structural block.
 */
struct StructuralBlockDefinition {
    std::string blockId;
    std::vector<std::string> siteIds;
};

enum class BlockOrientation {
    FORWARD,
    REVERSE,
    UNKNOWN
};

std::string blockOrientationName(BlockOrientation orientation);

/** HEAD precedes the first canonical site; TAIL follows the last one. */
enum class BlockEnd {
    HEAD,
    TAIL
};

std::string blockEndName(BlockEnd end);

struct OrientedBlockEndpoint {
    std::string blockId;
    BlockEnd end = BlockEnd::HEAD;

    bool operator<(const OrientedBlockEndpoint &other) const;
    bool operator==(const OrientedBlockEndpoint &other) const;
};

/** One exact projection of a structural block onto an extant occurrence path. */
struct ExtantBlockProjection {
    std::string projectionId;
    std::string taxon;
    std::string occurrencePath;
    std::string blockId;
    int64_t start0 = 0;
    int64_t end0 = 0;
    BlockOrientation orientation = BlockOrientation::UNKNOWN;
};

enum class BlockProjectionIssueKind {
    UNPROJECTED_BLOCK,
    UNKNOWN_ORIENTATION,
    DUPLICATE_BLOCK_COPY,
    OVERLAPPING_PROJECTIONS,
    SELF_ADJACENCY
};

std::string blockProjectionIssueKindName(BlockProjectionIssueKind kind);

struct BlockProjectionIssue {
    std::string issueId;
    BlockProjectionIssueKind kind = BlockProjectionIssueKind::UNPROJECTED_BLOCK;
    std::string taxon;
    std::string occurrencePath;
    std::vector<std::string> blockIds;
    std::vector<std::string> projectionIds;
    std::string detail;
};

/**
 * Read-only projections derived from the alignment graph or supplied by a
 * macro-synteny caller.  immutableCoreHash prevents accidental mixing with a
 * different or subsequently repaired graph.
 */
struct BlockProjectionSet {
    std::string immutableCoreHash;
    std::vector<ExtantBlockProjection> projections;
    std::vector<BlockProjectionIssue> issues;

    void validateAgainst(const AlignmentSiteGraph &graph) const;
};

struct BlockProjectionOptions {
    /** Require every site of a definition to have one resolved copy group. */
    bool requireCopyResolved = true;

    /** Require all sites in a definition to share that resolved copy group. */
    bool requireSingleCopyGroup = true;
};

struct AncestralAdjacencyConfig {
    std::string targetNode;
    BinaryParsimonyModel adjacencyModel;

    /** Include tied calls when auditing possible endpoint-degree conflicts. */
    bool ambiguousCallsCanConflict = true;
};

struct AncestralAdjacencyCall {
    std::string adjacencyId;
    OrientedBlockEndpoint endpoint1;
    OrientedBlockEndpoint endpoint2;
    BinaryInferenceResult presence;
    std::map<std::string, PresenceObservation> leafObservations;

    /** Same two blocks support more than one oriented endpoint pairing. */
    bool orientationConflict = false;

    /** This endpoint participates in another present or potentially present call. */
    bool endpointConflict = false;

    /** Present, non-ambiguous, and free of the two local conflict classes. */
    bool supportedWithoutLocalConflict = false;

    std::string decision;
};

/**
 * Candidate ancestral adjacencies only.  This object never contains a base
 * sequence and intentionally does not claim a chromosome ordering or a full
 * ancestral karyotype.  Telomere evidence, global matching, and macro-synteny
 * validation must be added before chromosome reconstruction is justified.
 */
struct AncestralAdjacencyReport {
    std::string immutableCoreHash;
    std::string targetNode;
    std::string inferenceScope = "CANDIDATE_BLOCK_ADJACENCIES_ONLY";
    std::vector<AncestralAdjacencyCall> calls;
    std::vector<BlockProjectionIssue> issues;

    void validateAgainst(const AlignmentSiteGraph &graph) const;
};

class AncestralAdjacencyError : public std::runtime_error {
public:
    explicit AncestralAdjacencyError(const std::string &message)
        : std::runtime_error(message) {}
};

class AncestralAdjacencyInference {
public:
    /**
     * Conservatively project exact, contiguous block definitions onto graph
     * paths.  A singleton or direction-symmetric definition is retained with
     * UNKNOWN orientation and cannot create an oriented adjacency.
     */
    static BlockProjectionSet projectBlocks(
        const AlignmentSiteGraph &graph,
        const std::vector<StructuralBlockDefinition> &definitions,
        const BlockProjectionOptions &options = BlockProjectionOptions());

    /**
     * Infer the presence of each observed oriented adjacency at targetNode.
     *
     * A leaf is PRESENT when it has one unambiguous occurrence of each block
     * and the requested ends are adjacent.  It is ABSENT when both unique
     * occurrences are callable but not adjacent, and MISSING when a block is
     * absent, duplicated, overlapping, or unoriented.  This conservative rule
     * keeps missing projection evidence distinct from rearrangement evidence.
     */
    static AncestralAdjacencyReport infer(
        const AlignmentSiteGraph &graph,
        const Phylogeny &tree,
        const BlockProjectionSet &projectionSet,
        const AncestralAdjacencyConfig &config);
};

class AncestralAdjacencyWriters {
public:
    static void writeProjections(const BlockProjectionSet &projectionSet,
                                 std::ostream &output);
    static void writeCalls(const AncestralAdjacencyReport &report,
                           std::ostream &output);
    static void writeIssues(const AncestralAdjacencyReport &report,
                            std::ostream &output);
};

}  // namespace trio
}  // namespace anchorwave
