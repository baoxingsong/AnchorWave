#include "Phylogeny.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace anchorwave {
namespace trio {

namespace {

Phylogeny::NodeId noNode() {
    return std::numeric_limits<Phylogeny::NodeId>::max();
}

bool nearlyEqual(double lhs, double rhs) {
    const double scale = std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
    return std::fabs(lhs - rhs) <= 1e-12 * scale;
}

template <std::size_t StateCount>
void validateCostMatrix(
    const std::array<std::array<double, StateCount>, StateCount>& costs,
    const char* modelName) {
    for (std::size_t parent = 0; parent < StateCount; ++parent) {
        for (std::size_t child = 0; child < StateCount; ++child) {
            const double cost = costs[parent][child];
            if (!std::isfinite(cost) || cost < 0.0) {
                std::ostringstream message;
                message << modelName << " transition cost [" << parent << "]["
                        << child << "] must be finite and non-negative";
                throw std::invalid_argument(message.str());
            }
        }
    }
}

template <std::size_t StateCount>
std::array<double, StateCount> sankoffAtNode(
    const Phylogeny& tree,
    Phylogeny::NodeId target,
    const std::vector<unsigned int>& allowedStateMasks,
    const std::array<std::array<double, StateCount>, StateCount>& transitionCost) {
    typedef std::array<double, StateCount> CostVector;

    std::function<CostVector(Phylogeny::NodeId, Phylogeny::NodeId)> visit;
    visit = [&](Phylogeny::NodeId node, Phylogeny::NodeId blocked) -> CostVector {
        CostVector result;
        for (std::size_t state = 0; state < StateCount; ++state) {
            const bool allowed =
                !tree.isLeaf(node) ||
                (allowedStateMasks[node] & (1u << static_cast<unsigned int>(state))) != 0u;
            result[state] = allowed ? 0.0 : std::numeric_limits<double>::infinity();
        }

        std::vector<Phylogeny::NodeId> neighbours;
        if (tree.hasParent(node)) {
            neighbours.push_back(tree.parent(node));
        }
        const std::vector<Phylogeny::NodeId>& childNodes = tree.children(node);
        neighbours.insert(neighbours.end(), childNodes.begin(), childNodes.end());

        for (std::size_t neighbourIndex = 0; neighbourIndex < neighbours.size();
             ++neighbourIndex) {
            const Phylogeny::NodeId neighbour = neighbours[neighbourIndex];
            if (neighbour == blocked) {
                continue;
            }

            const CostVector neighbourCosts = visit(neighbour, node);
            const bool nodeIsOriginalParent =
                tree.hasParent(neighbour) && tree.parent(neighbour) == node;

            for (std::size_t nodeState = 0; nodeState < StateCount; ++nodeState) {
                double bestNeighbourCost = std::numeric_limits<double>::infinity();
                for (std::size_t neighbourState = 0; neighbourState < StateCount;
                     ++neighbourState) {
                    const double edgeCost = nodeIsOriginalParent
                        ? transitionCost[nodeState][neighbourState]
                        : transitionCost[neighbourState][nodeState];
                    bestNeighbourCost = std::min(
                        bestNeighbourCost,
                        neighbourCosts[neighbourState] + edgeCost);
                }
                result[nodeState] += bestNeighbourCost;
            }
        }
        return result;
    };

    return visit(target, noNode());
}

template <std::size_t StateCount>
void summarizeCosts(const std::array<double, StateCount>& stateCosts,
                    std::vector<std::size_t>* optimalStateIndices,
                    double* bestCost,
                    double* secondBestCost,
                    double* confidenceMargin) {
    std::array<double, StateCount> sortedCosts = stateCosts;
    std::sort(sortedCosts.begin(), sortedCosts.end());

    *bestCost = sortedCosts[0];
    *secondBestCost = sortedCosts[1];
    optimalStateIndices->clear();
    for (std::size_t state = 0; state < StateCount; ++state) {
        if (nearlyEqual(stateCosts[state], *bestCost)) {
            optimalStateIndices->push_back(state);
        }
    }

    const double rawMargin = *secondBestCost - *bestCost;
    *confidenceMargin = nearlyEqual(*secondBestCost, *bestCost)
        ? 0.0
        : std::max(0.0, rawMargin);
}

unsigned int nucleotideMask(char symbol) {
    const char upper = static_cast<char>(
        std::toupper(static_cast<unsigned char>(symbol)));
    switch (upper) {
        case 'A': return 0x1u;
        case 'C': return 0x2u;
        case 'G': return 0x4u;
        case 'T':
        case 'U': return 0x8u;
        case 'R': return 0x1u | 0x4u;
        case 'Y': return 0x2u | 0x8u;
        case 'S': return 0x2u | 0x4u;
        case 'W': return 0x1u | 0x8u;
        case 'K': return 0x4u | 0x8u;
        case 'M': return 0x1u | 0x2u;
        case 'B': return 0x2u | 0x4u | 0x8u;
        case 'D': return 0x1u | 0x4u | 0x8u;
        case 'H': return 0x1u | 0x2u | 0x8u;
        case 'V': return 0x1u | 0x2u | 0x4u;
        case 'N':
        case 'X':
        case '?':
        case '.':
        case '-': return 0xFu;
        default: {
            std::ostringstream message;
            message << "invalid nucleotide observation symbol '" << symbol << "'";
            throw std::invalid_argument(message.str());
        }
    }
}

char iupacForMask(unsigned int mask) {
    switch (mask) {
        case 0x1u: return 'A';
        case 0x2u: return 'C';
        case 0x4u: return 'G';
        case 0x8u: return 'T';
        case 0x5u: return 'R';
        case 0xAu: return 'Y';
        case 0x6u: return 'S';
        case 0x9u: return 'W';
        case 0xCu: return 'K';
        case 0x3u: return 'M';
        case 0xEu: return 'B';
        case 0xDu: return 'D';
        case 0xBu: return 'H';
        case 0x7u: return 'V';
        case 0xFu: return 'N';
        default: throw std::logic_error("cannot encode an empty nucleotide state set");
    }
}

std::unordered_map<std::string, Phylogeny::NodeId> leafIndex(const Phylogeny& tree) {
    std::unordered_map<std::string, Phylogeny::NodeId> index;
    for (Phylogeny::NodeId node = 0; node < tree.nodeCount(); ++node) {
        if (tree.isLeaf(node)) {
            index.insert(std::make_pair(tree.nodeName(node), node));
        }
    }
    return index;
}

}  // namespace

class NewickParser {
public:
    explicit NewickParser(const std::string& input)
        : input_(input), position_(0) {}

    Phylogeny parse() {
        skipIgnored();
        if (atEnd()) {
            fail("tree is empty");
        }

        const Phylogeny::NodeId root = parseSubtree();
        skipIgnored();
        if (!atEnd() && input_[position_] == ';') {
            ++position_;
            skipIgnored();
        }
        if (!atEnd()) {
            fail("unexpected trailing content");
        }

        std::set<std::string> leaves;
        for (std::size_t node = 0; node < nodes_.size(); ++node) {
            if (nodes_[node].children.empty()) {
                if (nodes_[node].name.empty()) {
                    fail("every leaf must have a label");
                }
                if (!leaves.insert(nodes_[node].name).second) {
                    std::ostringstream message;
                    message << "duplicate leaf label '" << nodes_[node].name << "'";
                    fail(message.str());
                }
            }
        }

        return Phylogeny(std::move(nodes_), root);
    }

private:
    struct ParsedLength {
        bool present;
        double value;
    };

    const std::string& input_;
    std::size_t position_;
    std::vector<Phylogeny::Node> nodes_;

    bool atEnd() const {
        return position_ >= input_.size();
    }

    void fail(const std::string& reason) const {
        std::ostringstream message;
        message << "invalid Newick at position " << position_ << ": " << reason;
        throw std::invalid_argument(message.str());
    }

    void skipIgnored() {
        bool advanced = true;
        while (advanced) {
            advanced = false;
            while (!atEnd() &&
                   std::isspace(static_cast<unsigned char>(input_[position_]))) {
                ++position_;
                advanced = true;
            }
            if (!atEnd() && input_[position_] == '[') {
                advanced = true;
                int depth = 0;
                do {
                    const char current = input_[position_++];
                    if (current == '[') {
                        ++depth;
                    } else if (current == ']') {
                        --depth;
                    }
                    if (atEnd() && depth != 0) {
                        fail("unterminated comment");
                    }
                } while (depth != 0);
            }
        }
    }

    std::string parseOptionalLabel() {
        skipIgnored();
        if (atEnd()) {
            return std::string();
        }

        if (input_[position_] == '\'') {
            ++position_;
            std::string label;
            while (!atEnd()) {
                if (input_[position_] == '\'') {
                    if (position_ + 1 < input_.size() && input_[position_ + 1] == '\'') {
                        label.push_back('\'');
                        position_ += 2;
                        continue;
                    }
                    ++position_;
                    return label;
                }
                label.push_back(input_[position_++]);
            }
            fail("unterminated quoted label");
        }

        const char first = input_[position_];
        if (first == '(' || first == ')' || first == ',' || first == ':' ||
            first == ';' || first == '[' || first == ']') {
            return std::string();
        }

        const std::size_t begin = position_;
        while (!atEnd()) {
            const char current = input_[position_];
            if (std::isspace(static_cast<unsigned char>(current)) ||
                current == '(' || current == ')' || current == ',' || current == ':' ||
                current == ';' || current == '[' || current == ']') {
                break;
            }
            ++position_;
        }
        return input_.substr(begin, position_ - begin);
    }

    ParsedLength parseOptionalLength() {
        skipIgnored();
        if (atEnd() || input_[position_] != ':') {
            return ParsedLength{false, 0.0};
        }

        ++position_;
        skipIgnored();
        if (atEnd()) {
            fail("missing branch length after ':'");
        }

        errno = 0;
        const char* begin = input_.c_str() + position_;
        char* end = NULL;
        const double value = std::strtod(begin, &end);
        if (end == begin) {
            fail("invalid branch length");
        }
        position_ += static_cast<std::size_t>(end - begin);
        if (errno == ERANGE || !std::isfinite(value) || value < 0.0) {
            fail("branch length must be finite and non-negative");
        }
        return ParsedLength{true, value};
    }

    Phylogeny::NodeId appendNode(const std::string& name,
                                 const std::vector<Phylogeny::NodeId>& children,
                                 const ParsedLength& length) {
        Phylogeny::Node node;
        node.name = name;
        node.parent = noNode();
        node.children = children;
        node.hasBranchLength = length.present;
        node.branchLength = length.value;

        const Phylogeny::NodeId id = nodes_.size();
        nodes_.push_back(node);
        for (std::size_t i = 0; i < children.size(); ++i) {
            if (nodes_[children[i]].parent != noNode()) {
                fail("a node has more than one parent");
            }
            nodes_[children[i]].parent = id;
        }
        return id;
    }

    Phylogeny::NodeId parseSubtree() {
        skipIgnored();
        if (atEnd()) {
            fail("expected a subtree");
        }

        if (input_[position_] != '(') {
            const std::string leafName = parseOptionalLabel();
            if (leafName.empty()) {
                fail("expected a leaf label");
            }
            const ParsedLength length = parseOptionalLength();
            return appendNode(leafName, std::vector<Phylogeny::NodeId>(), length);
        }

        ++position_;
        skipIgnored();
        if (!atEnd() && input_[position_] == ')') {
            fail("internal node has no children");
        }

        std::vector<Phylogeny::NodeId> childNodes;
        childNodes.push_back(parseSubtree());
        while (true) {
            skipIgnored();
            if (atEnd()) {
                fail("unterminated internal node");
            }
            if (input_[position_] == ',') {
                ++position_;
                childNodes.push_back(parseSubtree());
                continue;
            }
            if (input_[position_] == ')') {
                ++position_;
                break;
            }
            fail("expected ',' or ')' after child subtree");
        }

        const std::string internalName = parseOptionalLabel();
        const ParsedLength length = parseOptionalLength();
        return appendNode(internalName, childNodes, length);
    }
};

NucleotideParsimonyModel::NucleotideParsimonyModel() {
    for (std::size_t parent = 0; parent < 4; ++parent) {
        for (std::size_t child = 0; child < 4; ++child) {
            transitionCost[parent][child] = parent == child ? 0.0 : 1.0;
        }
    }
}

NucleotideParsimonyModel NucleotideParsimonyModel::equalCost(
    double mismatchCost,
    double matchCost) {
    NucleotideParsimonyModel model;
    for (std::size_t parent = 0; parent < 4; ++parent) {
        for (std::size_t child = 0; child < 4; ++child) {
            model.transitionCost[parent][child] =
                parent == child ? matchCost : mismatchCost;
        }
    }
    return model;
}

BinaryParsimonyModel::BinaryParsimonyModel() {
    transitionCost[0][0] = 0.0;
    transitionCost[0][1] = 1.0;
    transitionCost[1][0] = 1.0;
    transitionCost[1][1] = 0.0;
}

BinaryParsimonyModel BinaryParsimonyModel::gainLoss(
    double gainCost,
    double lossCost,
    double retainAbsentCost,
    double retainPresentCost) {
    BinaryParsimonyModel model;
    model.transitionCost[0][0] = retainAbsentCost;
    model.transitionCost[0][1] = gainCost;
    model.transitionCost[1][0] = lossCost;
    model.transitionCost[1][1] = retainPresentCost;
    return model;
}

bool NucleotideInferenceResult::isAmbiguous() const {
    return optimalStates.size() != 1u;
}

bool BinaryInferenceResult::isAmbiguous() const {
    return optimalStates.size() != 1u;
}

Phylogeny::Phylogeny(std::vector<Node> nodes, NodeId root)
    : nodes_(std::move(nodes)), root_(root) {}

Phylogeny Phylogeny::fromNewick(const std::string& newick) {
    return NewickParser(newick).parse();
}

Phylogeny::NodeId Phylogeny::root() const {
    return root_;
}

std::size_t Phylogeny::nodeCount() const {
    return nodes_.size();
}

std::size_t Phylogeny::leafCount() const {
    std::size_t count = 0;
    for (std::size_t node = 0; node < nodes_.size(); ++node) {
        if (nodes_[node].children.empty()) {
            ++count;
        }
    }
    return count;
}

void Phylogeny::validateNodeId(NodeId node) const {
    if (node >= nodes_.size()) {
        throw std::out_of_range("phylogeny node ID is out of range");
    }
}

void Phylogeny::validateInternalNode(NodeId node) const {
    validateNodeId(node);
    if (isLeaf(node)) {
        throw std::invalid_argument("ancestral inference target must be an internal node");
    }
}

bool Phylogeny::isLeaf(NodeId node) const {
    validateNodeId(node);
    return nodes_[node].children.empty();
}

bool Phylogeny::isInternal(NodeId node) const {
    return !isLeaf(node);
}

bool Phylogeny::hasParent(NodeId node) const {
    validateNodeId(node);
    return nodes_[node].parent != noNode();
}

Phylogeny::NodeId Phylogeny::parent(NodeId node) const {
    validateNodeId(node);
    if (!hasParent(node)) {
        throw std::logic_error("the root node has no parent");
    }
    return nodes_[node].parent;
}

const std::vector<Phylogeny::NodeId>& Phylogeny::children(NodeId node) const {
    validateNodeId(node);
    return nodes_[node].children;
}

const std::string& Phylogeny::nodeName(NodeId node) const {
    validateNodeId(node);
    return nodes_[node].name;
}

bool Phylogeny::hasBranchLength(NodeId node) const {
    validateNodeId(node);
    return nodes_[node].hasBranchLength;
}

double Phylogeny::branchLength(NodeId node) const {
    validateNodeId(node);
    if (!nodes_[node].hasBranchLength) {
        throw std::logic_error("node has no branch length");
    }
    return nodes_[node].branchLength;
}

std::vector<std::string> Phylogeny::leafNames() const {
    std::vector<std::string> names;
    names.reserve(leafCount());
    for (std::size_t node = 0; node < nodes_.size(); ++node) {
        if (nodes_[node].children.empty()) {
            names.push_back(nodes_[node].name);
        }
    }
    return names;
}

Phylogeny::NodeId Phylogeny::internalNodeId(const std::string& name) const {
    if (name.empty()) {
        throw std::invalid_argument("internal node name must not be empty");
    }

    NodeId found = noNode();
    bool leafHasName = false;
    for (NodeId node = 0; node < nodes_.size(); ++node) {
        if (nodes_[node].name != name) {
            continue;
        }
        if (nodes_[node].children.empty()) {
            leafHasName = true;
            continue;
        }
        if (found != noNode()) {
            throw std::invalid_argument("internal node name is ambiguous: " + name);
        }
        found = node;
    }

    if (found == noNode()) {
        if (leafHasName) {
            throw std::invalid_argument("name identifies a leaf, not an internal node: " + name);
        }
        throw std::invalid_argument("internal node name was not found: " + name);
    }
    return found;
}

NucleotideInferenceResult Phylogeny::inferNucleotide(
    const std::string& targetInternalName,
    const std::map<std::string, char>& leafObservations,
    const NucleotideParsimonyModel& model) const {
    return inferNucleotide(internalNodeId(targetInternalName), leafObservations, model);
}

NucleotideInferenceResult Phylogeny::inferNucleotide(
    NodeId targetInternalNode,
    const std::map<std::string, char>& leafObservations,
    const NucleotideParsimonyModel& model) const {
    validateInternalNode(targetInternalNode);
    validateCostMatrix(model.transitionCost, "nucleotide model");

    std::vector<unsigned int> masks(nodeCount(), 0xFu);
    const std::unordered_map<std::string, NodeId> leaves = leafIndex(*this);
    for (std::map<std::string, char>::const_iterator observation =
             leafObservations.begin();
         observation != leafObservations.end(); ++observation) {
        const std::unordered_map<std::string, NodeId>::const_iterator leaf =
            leaves.find(observation->first);
        if (leaf == leaves.end()) {
            throw std::invalid_argument(
                "nucleotide observation names an unknown leaf: " + observation->first);
        }
        masks[leaf->second] = nucleotideMask(observation->second);
    }

    NucleotideInferenceResult result;
    result.stateCosts = sankoffAtNode<4>(
        *this, targetInternalNode, masks, model.transitionCost);

    std::vector<std::size_t> optimalIndices;
    summarizeCosts(result.stateCosts,
                   &optimalIndices,
                   &result.bestCost,
                   &result.secondBestCost,
                   &result.confidenceMargin);

    static const char states[4] = {'A', 'C', 'G', 'T'};
    unsigned int optimalMask = 0u;
    for (std::size_t i = 0; i < optimalIndices.size(); ++i) {
        result.optimalStates.push_back(states[optimalIndices[i]]);
        optimalMask |= 1u << static_cast<unsigned int>(optimalIndices[i]);
    }
    result.selectedState = result.optimalStates.front();
    result.ambiguityCode = iupacForMask(optimalMask);
    return result;
}

BinaryInferenceResult Phylogeny::inferPresence(
    const std::string& targetInternalName,
    const std::map<std::string, PresenceObservation>& leafObservations,
    const BinaryParsimonyModel& model) const {
    return inferPresence(internalNodeId(targetInternalName), leafObservations, model);
}

BinaryInferenceResult Phylogeny::inferPresence(
    NodeId targetInternalNode,
    const std::map<std::string, PresenceObservation>& leafObservations,
    const BinaryParsimonyModel& model) const {
    validateInternalNode(targetInternalNode);
    validateCostMatrix(model.transitionCost, "binary model");

    std::vector<unsigned int> masks(nodeCount(), 0x3u);
    const std::unordered_map<std::string, NodeId> leaves = leafIndex(*this);
    for (std::map<std::string, PresenceObservation>::const_iterator observation =
             leafObservations.begin();
         observation != leafObservations.end(); ++observation) {
        const std::unordered_map<std::string, NodeId>::const_iterator leaf =
            leaves.find(observation->first);
        if (leaf == leaves.end()) {
            throw std::invalid_argument(
                "presence observation names an unknown leaf: " + observation->first);
        }

        switch (observation->second) {
            case PresenceObservation::ABSENT:
                masks[leaf->second] = 0x1u;
                break;
            case PresenceObservation::PRESENT:
                masks[leaf->second] = 0x2u;
                break;
            case PresenceObservation::MISSING:
                masks[leaf->second] = 0x3u;
                break;
            default:
                throw std::invalid_argument("invalid presence observation value");
        }
    }

    BinaryInferenceResult result;
    result.stateCosts = sankoffAtNode<2>(
        *this, targetInternalNode, masks, model.transitionCost);

    std::vector<std::size_t> optimalIndices;
    summarizeCosts(result.stateCosts,
                   &optimalIndices,
                   &result.bestCost,
                   &result.secondBestCost,
                   &result.confidenceMargin);
    for (std::size_t i = 0; i < optimalIndices.size(); ++i) {
        result.optimalStates.push_back(static_cast<int>(optimalIndices[i]));
    }
    result.selectedState = result.optimalStates.front();
    return result;
}

}  // namespace trio
}  // namespace anchorwave
