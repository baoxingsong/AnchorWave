#ifndef ANCHORWAVE_TRIO_PHYLOGENY_H
#define ANCHORWAVE_TRIO_PHYLOGENY_H

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

/**
 * A substitution-cost model used only for nucleotide ancestral inference.
 *
 * States are ordered A, C, G, T.  transitionCost[parent][child] is the
 * parsimony cost assigned to one directed tree edge.  This model is
 * deliberately independent of all sequence-alignment scores and gap costs.
 */
struct NucleotideParsimonyModel {
    std::array<std::array<double, 4>, 4> transitionCost;

    /** Construct the usual model with zero matches and unit mismatches. */
    NucleotideParsimonyModel();

    /** Construct a model with a common match and mismatch cost. */
    static NucleotideParsimonyModel equalCost(double mismatchCost = 1.0,
                                               double matchCost = 0.0);
};

/**
 * A gain/loss model used only for run-level presence/absence inference.
 *
 * States are ordered absent (0), present (1).  transitionCost[parent][child]
 * is directed, so gain and loss can have different costs.  These costs are
 * not derived from affine alignment penalties.
 */
struct BinaryParsimonyModel {
    std::array<std::array<double, 2>, 2> transitionCost;

    /** Construct a symmetric unit gain/loss model. */
    BinaryParsimonyModel();

    /** Construct a model with independently specified gain and loss costs. */
    static BinaryParsimonyModel gainLoss(double gainCost,
                                         double lossCost,
                                         double retainAbsentCost = 0.0,
                                         double retainPresentCost = 0.0);
};

/** Observed state for a homologous sequence run at one leaf. */
enum class PresenceObservation {
    ABSENT = 0,
    PRESENT = 1,
    MISSING = 2
};

/** Result of nucleotide Sankoff inference at one internal node. */
struct NucleotideInferenceResult {
    /** Global minimum tree cost when the target is fixed to A, C, G, or T. */
    std::array<double, 4> stateCosts;

    /** All minimum-cost states, always ordered A, C, G, T. */
    std::vector<char> optimalStates;

    /** First optimal state in A, C, G, T order, for deterministic tie handling. */
    char selectedState;

    /** IUPAC encoding of the complete optimal-state set. */
    char ambiguityCode;

    double bestCost;

    /** Second-lowest candidate-state cost, including an equal-cost tie. */
    double secondBestCost;

    /** secondBestCost - bestCost; therefore zero whenever the optimum is tied. */
    double confidenceMargin;

    bool isAmbiguous() const;
};

/** Result of binary Sankoff inference at one internal node. */
struct BinaryInferenceResult {
    /** Global minimum tree cost when the target is fixed to absent or present. */
    std::array<double, 2> stateCosts;

    /** All minimum-cost states, always ordered absent (0), present (1). */
    std::vector<int> optimalStates;

    /** First optimal state; absent wins an exact tie only as a deterministic API value. */
    int selectedState;

    double bestCost;
    double secondBestCost;
    double confidenceMargin;

    bool isAmbiguous() const;
};

/**
 * Immutable rooted phylogeny parsed from Newick.
 *
 * The parser accepts named or unnamed internal nodes, quoted labels, comments,
 * whitespace, and optional non-negative branch lengths.  Leaf labels must be
 * present and unique.  Branch lengths are retained for later likelihood-model
 * use but are not silently applied to Sankoff costs.
 *
 * Sankoff inference may target any internal node, not only the root.  It uses
 * evidence from every observed leaf on both sides of the target and therefore
 * supports any number of outgroups.  Missing leaf observations may either be
 * omitted from the observation map or represented explicitly.
 */
class Phylogeny {
public:
    typedef std::size_t NodeId;

    /** Parse one Newick tree.  A trailing semicolon is optional. */
    static Phylogeny fromNewick(const std::string& newick);

    NodeId root() const;
    std::size_t nodeCount() const;
    std::size_t leafCount() const;

    bool isLeaf(NodeId node) const;
    bool isInternal(NodeId node) const;
    bool hasParent(NodeId node) const;
    NodeId parent(NodeId node) const;
    const std::vector<NodeId>& children(NodeId node) const;

    /** Return the parsed label, or an empty string for an unnamed internal node. */
    const std::string& nodeName(NodeId node) const;

    bool hasBranchLength(NodeId node) const;
    double branchLength(NodeId node) const;

    /** Leaf labels in deterministic parse order. */
    std::vector<std::string> leafNames() const;

    /**
     * Find a uniquely named internal node.
     *
     * Throws std::invalid_argument when the name is absent, names a leaf, or
     * is shared by more than one internal node.
     */
    NodeId internalNodeId(const std::string& name) const;

    /**
     * Infer an ancestral nucleotide at a named internal node.
     *
     * Observation symbols are case-insensitive IUPAC nucleotide codes.
     * '?', '.', '-', 'N', and 'X' provide no nucleotide-state constraint.
     * A gap's biological presence/absence must be inferred separately with
     * inferPresence(); '-' here means only that no nucleotide was observed.
     * Leaves omitted from the map are treated as missing.  Unknown taxa and
     * invalid symbols are rejected.
     */
    NucleotideInferenceResult inferNucleotide(
        const std::string& targetInternalName,
        const std::map<std::string, char>& leafObservations,
        const NucleotideParsimonyModel& model = NucleotideParsimonyModel()) const;

    /** Same as the named overload, with a previously validated internal node ID. */
    NucleotideInferenceResult inferNucleotide(
        NodeId targetInternalNode,
        const std::map<std::string, char>& leafObservations,
        const NucleotideParsimonyModel& model = NucleotideParsimonyModel()) const;

    /**
     * Infer ancestral presence/absence with an independently parameterized
     * binary gain/loss model.  Omitted leaves are treated as MISSING.
     */
    BinaryInferenceResult inferPresence(
        const std::string& targetInternalName,
        const std::map<std::string, PresenceObservation>& leafObservations,
        const BinaryParsimonyModel& model = BinaryParsimonyModel()) const;

    /** Same as the named overload, with a previously validated internal node ID. */
    BinaryInferenceResult inferPresence(
        NodeId targetInternalNode,
        const std::map<std::string, PresenceObservation>& leafObservations,
        const BinaryParsimonyModel& model = BinaryParsimonyModel()) const;

private:
    struct Node {
        std::string name;
        NodeId parent;
        std::vector<NodeId> children;
        bool hasBranchLength;
        double branchLength;
    };

    std::vector<Node> nodes_;
    NodeId root_;

    explicit Phylogeny(std::vector<Node> nodes, NodeId root);

    void validateNodeId(NodeId node) const;
    void validateInternalNode(NodeId node) const;

    friend class NewickParser;
};

}  // namespace trio
}  // namespace anchorwave

#endif  // ANCHORWAVE_TRIO_PHYLOGENY_H
