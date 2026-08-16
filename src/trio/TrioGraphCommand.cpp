#include "src/trio/TrioGraphCommand.h"

#include "src/trio/Phylogeny.h"
#include "src/trio/evolution/AncestralAdjacency.h"
#include "src/trio/evolution/AncestralInference.h"
#include "src/trio/impl/CopyRelationshipResolver.h"
#include "src/trio/impl/GraphRepairEngine.h"
#include "src/trio/io/BlockProjectionReader.h"
#include "src/trio/io/CopyConstraintReader.h"
#include "src/trio/io/FastaEvidenceValidator.h"
#include "src/trio/io/GraphWriters.h"
#include "src/trio/io/PairwiseManifestReader.h"
#include "src/trio/io/TaxonManifestReader.h"
#include "src/trio/model/AlignmentEvidence.h"
#include "src/trio/model/AlignmentSiteGraph.h"

#include <cstdio>
#include <fstream>
#include <functional>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace anchorwave {
namespace trio {
namespace {

struct CommandOptions {
    TrioConfig config;
    std::string ancestorNode;
    std::string blockProjectionsPath;
    bool validateOnly = false;
    bool skipAncestry = false;
    bool includeAmbiguousPresence = false;
    bool strictCopyMode = false;
    bool localRepair = true;
    bool validateFastaEvidence = true;
    GraphRepairOptions repairOptions;
};

std::string usageText() {
    return
        "Usage: anchorwave trioGraphAli --taxa taxa.tsv --pairwise-manifest pairwise.tsv\n"
        "       --species-tree species.nwk --output-prefix RUN [options]\n\n"
        "Required:\n"
        "  --taxa FILE                 versioned taxon manifest\n"
        "  --pairwise-manifest FILE    all three primary-trio MAF edges; 2 per extra outgroup\n"
        "  --species-tree FILE         rooted Newick tree containing every manifest taxon\n"
        "  --output-prefix PREFIX      output prefix; parent directory must exist\n\n"
        "Copy resolution:\n"
        "  --copy-relations FILE       versioned copy GROUP/MEMBER/COUNT TSV (optional)\n"
        "  --copy-number SPEC          repeatable TAXON=N or FAMILY:TAXON=N capacity\n"
        "  --copy-mode MODE            constrained (default) or strict\n\n"
        "Evidence and inference:\n"
        "  --pairwise-scope MODE       triangles (default) or complete\n"
        "  --ancestor-node NAME        infer this named internal node\n"
        "  --block-projections FILE    independent extant macro-block projections\n"
        "  --include-ambiguous-presence emit ambiguous-presence sites in ancestor blocks\n"
        "  --no-local-repair           classify discordance but do not run bounded POA\n"
        "  --skip-fasta-validation     do not compare MAF spelling/lengths with FASTA\n"
        "  --poa-min-delta FLOAT       minimum alignment-only repair improvement (default 0)\n"
        "  --poa-max-core-sites N      maximum mutable graph sites (default 10000)\n"
        "  --poa-max-dp-cells N        maximum cells in one local DP (default 5000000)\n"
        "  --poa-match FLOAT           alignment match score (default 2)\n"
        "  --poa-mismatch FLOAT        alignment mismatch score (default -4)\n"
        "  --poa-open1/--poa-extend1 FLOAT  first affine gap scores\n"
        "  --poa-open2/--poa-extend2 FLOAT  second affine gap scores\n"
        "  --skip-ancestry             build only the immutable extant graph\n"
        "  --validate-input-paths      check all manifest paths before parsing MAFs\n"
        "  --validate-only             validate manifests/tree/MAFs/copies; write no outputs\n"
        "  -h, --help                  show this help\n\n"
        "The GFA plus graph metadata are authoritative. MAF files are derived projections.\n";
}

std::string requireValue(int &index, int argc, char **argv,
                         const std::string &option) {
    if (index + 1 >= argc) throw std::invalid_argument(option + " requires a value");
    return argv[++index];
}

double parseDoubleValue(const std::string &text, const std::string &option) {
    std::size_t consumed = 0;
    double value = 0.0;
    try { value = std::stod(text, &consumed); }
    catch (const std::exception &) {
        throw std::invalid_argument(option + " requires a finite number");
    }
    if (consumed != text.size() || !std::isfinite(value))
        throw std::invalid_argument(option + " requires a finite number");
    return value;
}

std::size_t parseSizeValue(const std::string &text, const std::string &option) {
    if (text.empty() || text[0] == '-') {
        throw std::invalid_argument(option + " requires a non-negative integer");
    }
    std::size_t consumed = 0;
    unsigned long long value = 0;
    try { value = std::stoull(text, &consumed); }
    catch (const std::exception &) {
        throw std::invalid_argument(option + " requires a non-negative integer");
    }
    if (consumed != text.size() ||
        value > static_cast<unsigned long long>(
                    std::numeric_limits<std::size_t>::max()))
        throw std::invalid_argument(option + " requires a non-negative integer");
    return static_cast<std::size_t>(value);
}

CommandOptions parseOptions(int argc, char **argv) {
    CommandOptions result;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "-h" || option == "--help") {
            std::cout << usageText();
            throw std::runtime_error("__help__");
        } else if (option == "--taxa") {
            result.config.taxaManifestPath = requireValue(i, argc, argv, option);
        } else if (option == "--pairwise-manifest") {
            result.config.pairwiseManifestPath = requireValue(i, argc, argv, option);
        } else if (option == "--species-tree") {
            result.config.speciesTreePath = requireValue(i, argc, argv, option);
        } else if (option == "--copy-relations") {
            result.config.copyRelationsPath = requireValue(i, argc, argv, option);
        } else if (option == "--copy-number") {
            result.config.defaultCopyCapacitySpecs.push_back(
                requireValue(i, argc, argv, option));
        } else if (option == "--copy-mode") {
            const std::string mode = requireValue(i, argc, argv, option);
            if (mode == "strict") result.strictCopyMode = true;
            else if (mode != "constrained")
                throw std::invalid_argument("--copy-mode must be constrained or strict");
        } else if (option == "--pairwise-scope") {
            const std::string scope = requireValue(i, argc, argv, option);
            if (scope == "triangles") result.config.pairwiseScope = PairwiseScope::TRIANGLES;
            else if (scope == "complete") result.config.pairwiseScope = PairwiseScope::COMPLETE;
            else throw std::invalid_argument("--pairwise-scope must be triangles or complete");
        } else if (option == "--output-prefix") {
            result.config.outputPrefix = requireValue(i, argc, argv, option);
        } else if (option == "--ancestor-node") {
            result.ancestorNode = requireValue(i, argc, argv, option);
        } else if (option == "--block-projections") {
            result.blockProjectionsPath = requireValue(i, argc, argv, option);
        } else if (option == "--validate-input-paths") {
            result.config.validateInputPaths = true;
        } else if (option == "--validate-only") {
            result.validateOnly = true;
        } else if (option == "--skip-ancestry") {
            result.skipAncestry = true;
        } else if (option == "--include-ambiguous-presence") {
            result.includeAmbiguousPresence = true;
        } else if (option == "--no-local-repair") {
            result.localRepair = false;
        } else if (option == "--skip-fasta-validation") {
            result.validateFastaEvidence = false;
        } else if (option == "--poa-min-delta") {
            result.repairOptions.poa.minimumScoreDelta = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-max-core-sites") {
            result.repairOptions.maximumCoreSites = parseSizeValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-max-dp-cells") {
            result.repairOptions.poa.maxPairwiseDpCells = parseSizeValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-match") {
            result.repairOptions.poa.scoring.match = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-mismatch") {
            result.repairOptions.poa.scoring.mismatch = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-open1") {
            result.repairOptions.poa.scoring.gapOpen1 = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-extend1") {
            result.repairOptions.poa.scoring.gapExtend1 = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-open2") {
            result.repairOptions.poa.scoring.gapOpen2 = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else if (option == "--poa-extend2") {
            result.repairOptions.poa.scoring.gapExtend2 = parseDoubleValue(
                requireValue(i, argc, argv, option), option);
        } else {
            throw std::invalid_argument("unknown trioGraphAli option: " + option);
        }
    }
    if (result.config.taxaManifestPath.empty() ||
        result.config.pairwiseManifestPath.empty() ||
        result.config.speciesTreePath.empty()) {
        throw std::invalid_argument(
            "--taxa, --pairwise-manifest, and --species-tree are required");
    }
    if (!result.validateOnly && result.config.outputPrefix.empty()) {
        throw std::invalid_argument("--output-prefix is required unless --validate-only is used");
    }
    if (!result.skipAncestry && !result.validateOnly && result.ancestorNode.empty()) {
        throw std::invalid_argument(
            "--ancestor-node is required for ancestry output (or use --skip-ancestry)");
    }
    if (result.skipAncestry && !result.blockProjectionsPath.empty()) {
        throw std::invalid_argument(
            "--block-projections cannot be combined with --skip-ancestry");
    }
    return result;
}

std::string readTextFile(const std::string &path) {
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error(path + ": cannot open file");
    std::ostringstream value;
    value << input.rdbuf();
    if (!input.eof() && input.fail()) throw std::runtime_error(path + ": read failure");
    return value.str();
}

std::string resolveManifestEntry(const std::string &manifestPath,
                                 const std::string &entry) {
    if (entry.empty() || entry == "." || entry[0] == '/') return entry;
    const std::size_t separator = manifestPath.find_last_of("/\\");
    if (separator == std::string::npos) return entry;
    return manifestPath.substr(0, separator + 1) + entry;
}

void validateTreeTaxa(const Phylogeny &tree, const TaxonManifest &taxa,
                      const std::string &treePath) {
    const std::vector<std::string> leaves = tree.leafNames();
    const std::set<std::string> leafSet(leaves.begin(), leaves.end());
    const std::vector<std::string> declared = taxa.allTaxonIds();
    std::vector<std::string> missing;
    for (const std::string &taxon : declared) {
        if (!leafSet.count(taxon)) missing.push_back(taxon);
    }
    if (!missing.empty()) {
        std::ostringstream message;
        message << treePath << ": species tree is missing manifest taxa:";
        for (const std::string &taxon : missing) message << ' ' << taxon;
        throw std::runtime_error(message.str());
    }
}

template <typename Writer>
void writeAtomic(const std::string &path, Writer writer) {
    const std::string temporary = path + ".tmp";
    {
        std::ofstream output(temporary.c_str(), std::ios::binary);
        if (!output) throw std::runtime_error(path + ": cannot create temporary output");
        writer(output);
        output.close();
        if (!output) {
            std::remove(temporary.c_str());
            throw std::runtime_error(path + ": output write failed");
        }
    }
    if (std::rename(temporary.c_str(), path.c_str()) != 0) {
        std::remove(temporary.c_str());
        throw std::runtime_error(path + ": cannot atomically install output");
    }
}

void writeCopyAudit(const CopyResolutionResult &copies, std::ostream &output) {
    output << "record_type\tresidue_or_evidence\tselected\tcopy_group"
              "\tclassification\treason_or_provenance\n";
    for (const auto &entry : copies.assignments) {
        output << "ASSIGNMENT\t" << entry.first.canonicalString() << "\t1\t"
               << entry.second.copyGroup << "\t"
               << (entry.second.hard ? "HARD" : "INFERRED_OR_SOFT") << '\t'
               << entry.second.provenance << '\n';
    }
    for (const CopyResolutionAudit &audit : copies.audit) {
        output << "EVIDENCE\t" << audit.evidenceId << '\t'
               << (audit.selected ? 1 : 0) << "\t.\t" << audit.classification
               << '\t' << audit.reason << '\n';
    }
    for (const ResidueCopyAlternative &alternative : copies.alternatives) {
        output << "ALTERNATIVE\t" << alternative.residue.canonicalString()
               << "\t0\t";
        for (std::size_t i = 0; i < alternative.copyGroups.size(); ++i) {
            if (i) output << ',';
            output << alternative.copyGroups[i] << ':' << alternative.scores[i];
        }
        output << "\tSOFT_TIE\tretained_not_forced\n";
    }
}

void writeDiscordance(const AlignmentSiteGraph &graph, std::ostream &output) {
    output << "conflict_id\tclass\tpoa_eligible\tdisposition\tevidence_ids\tresidues\n";
    for (const GraphConflict &conflict : graph.conflicts()) {
        output << conflict.conflictId << '\t'
               << consistencyClassName(conflict.classification) << '\t'
               << (conflict.poaEligible ? 1 : 0) << '\t' << conflict.disposition << '\t';
        for (std::size_t i = 0; i < conflict.evidenceIds.size(); ++i) {
            if (i) output << ',';
            output << conflict.evidenceIds[i];
        }
        output << '\t';
        for (std::size_t i = 0; i < conflict.residues.size(); ++i) {
            if (i) output << ',';
            output << conflict.residues[i].canonicalString();
        }
        output << '\n';
    }
}

void writeRepairAudit(const std::vector<GraphRepairAudit> &audits,
                      std::ostream &output) {
    output << "conflict_id\trepair_id\tdisposition\tleft_site\tright_site"
              "\told_core_sites\tnew_core_sites\toutside_hash_before"
              "\toutside_hash_after\tbaseline_score\trepaired_score\tscore_delta\n";
    for (const GraphRepairAudit &audit : audits) {
        output << audit.conflictId << '\t'
               << (audit.repairId.empty() ? "." : audit.repairId) << '\t'
               << audit.disposition << '\t'
               << (audit.immutableLeftSite.empty() ? "." : audit.immutableLeftSite) << '\t'
               << (audit.immutableRightSite.empty() ? "." : audit.immutableRightSite) << '\t';
        for (std::size_t i = 0; i < audit.oldCoreSites.size(); ++i) {
            if (i) output << ',';
            output << audit.oldCoreSites[i];
        }
        output << '\t';
        for (std::size_t i = 0; i < audit.newCoreSites.size(); ++i) {
            if (i) output << ',';
            output << audit.newCoreSites[i];
        }
        output << '\t' << (audit.outsideHashBefore.empty() ? "." : audit.outsideHashBefore)
               << '\t' << (audit.outsideHashAfter.empty() ? "." : audit.outsideHashAfter)
               << '\t' << audit.baselineScore << '\t' << audit.repairedScore
               << '\t' << audit.scoreDelta << '\n';
    }
}

std::string copyStatusName(CopyResolutionStatus status) {
    switch (status) {
        case CopyResolutionStatus::RESOLVED: return "RESOLVED";
        case CopyResolutionStatus::ALTERNATIVES: return "ALTERNATIVES";
        case CopyResolutionStatus::UNRESOLVED: return "UNRESOLVED";
    }
    return "UNKNOWN";
}

int run(const CommandOptions &options) {
    ManifestReaderOptions taxonOptions;
    taxonOptions.validatePaths = options.config.validateInputPaths;
    const TaxonManifest taxa = TaxonManifestReader(taxonOptions).read(
        options.config.taxaManifestPath);
    PairwiseManifestReaderOptions pairOptions;
    pairOptions.validatePaths = options.config.validateInputPaths;
    pairOptions.scope = options.config.pairwiseScope;
    const PairwiseManifest pairs = PairwiseManifestReader(pairOptions).read(
        options.config.pairwiseManifestPath, taxa);
    const CopyConstraintSet constraints = CopyConstraintReader(taxonOptions).read(
        options.config.copyRelationsPath, taxa,
        options.config.defaultCopyCapacitySpecs);
    const Phylogeny tree = Phylogeny::fromNewick(
        readTextFile(options.config.speciesTreePath));
    validateTreeTaxa(tree, taxa, options.config.speciesTreePath);
    if (!options.ancestorNode.empty()) tree.internalNodeId(options.ancestorNode);

    std::vector<PairwiseMafInput> mafInputs;
    for (const PairwiseManifestRecord &record : pairs.records) {
        PairwiseMafInput input;
        input.leftTaxon = record.sourceTaxonA;
        input.rightTaxon = record.sourceTaxonB;
        input.mafPath = resolveManifestEntry(options.config.pairwiseManifestPath,
                                             record.maf);
        input.runId = record.pair.stableId();
        input.sourceWeight = record.weight;
        mafInputs.push_back(input);
    }
    const AlignmentEvidenceSet evidence = PairwiseEvidenceNormalizer::normalize(mafInputs);
    FastaEvidenceValidationCounts fastaValidation;
    if (options.validateFastaEvidence) {
        FastaEvidenceValidator::TaxonFastaMap fastaByTaxon;
        for (const TaxonRecord &record : taxa.records) {
            fastaByTaxon[record.taxonId] = record.fasta;
        }
        fastaValidation = FastaEvidenceValidator::validate(evidence, fastaByTaxon);
    }
    CopyResolverOptions resolverOptions;
    resolverOptions.strictExplicitGroups = options.strictCopyMode;
    const CopyResolutionResult copies = CopyRelationshipResolver::resolve(
        evidence, constraints, resolverOptions);

    AlignmentGraphBuildOptions graphOptions;
    const std::vector<std::string> ingroups = taxa.ingroupIds();
    graphOptions.ingroup1 = ingroups[0];
    graphOptions.ingroup2 = ingroups[1];
    graphOptions.primaryOutgroup = taxa.primaryOutgroupId();
    graphOptions.copyAssignments = copies.assignments;
    graphOptions.allowedEvidenceIds = copies.selectedEvidenceIds;
    graphOptions.restrictToAllowedEvidence = true;
    graphOptions.requireResolvedCopyAssignments = true;
    AlignmentSiteGraph graph = AlignmentSiteGraphBuilder::build(evidence, graphOptions);
    std::vector<GraphRepairAudit> repairAudits;
    if (options.localRepair) {
        GraphRepairResult repaired = GraphRepairEngine::repairEligibleConflicts(
            graph, options.repairOptions);
        graph = repaired.graph;
        repairAudits = repaired.audit;
    }
    BlockProjectionSet structuralProjections;
    if (!options.blockProjectionsPath.empty()) {
        structuralProjections = BlockProjectionReader::read(
            options.blockProjectionsPath, taxa, graph);
    }
    if (options.validateOnly) {
        std::cout << "trioGraphAli validation passed: taxa=" << taxa.records.size()
                  << " pairs=" << pairs.records.size()
                  << " residues=" << evidence.residues.size()
                  << " fasta_residues_validated="
                  << fastaValidation.residueObservationsValidated
                  << " sites=" << graph.sites().size()
                  << " conflicts=" << graph.conflicts().size()
                  << " structural_projections="
                  << structuralProjections.projections.size() << '\n';
        return 0;
    }

    const std::string prefix = options.config.outputPrefix;
    writeAtomic(prefix + ".graph.gfa",
                [&](std::ostream &out) { GfaWriter::write(graph, out); });
    writeAtomic(prefix + ".graph.meta.tsv",
                [&](std::ostream &out) { GfaWriter::writeMetadata(graph, out); });
    MafExportResult mafResult;
    writeAtomic(prefix + ".extant.maf", [&](std::ostream &out) {
        mafResult = MafGraphExporter::write(graph, out);
    });
    writeAtomic(prefix + ".extant.maf.omissions.tsv", [&](std::ostream &out) {
        MafGraphExporter::writeOmissions(mafResult, graph.coreHash(), out);
    });
    writeAtomic(prefix + ".copy-resolution.tsv",
                [&](std::ostream &out) { writeCopyAudit(copies, out); });
    writeAtomic(prefix + ".discordance.tsv",
                [&](std::ostream &out) { writeDiscordance(graph, out); });
    writeAtomic(prefix + ".repair-audit.tsv",
                [&](std::ostream &out) { writeRepairAudit(repairAudits, out); });

    std::size_t ancestorBlocks = 0;
    std::size_t childBlocks = 0;
    std::size_t adjacencyCalls = 0;
    std::size_t supportedAdjacencies = 0;
    std::size_t adjacencyIssues = 0;
    if (!options.skipAncestry) {
        EvolutionInferenceConfig evolution;
        evolution.targetNode = options.ancestorNode;
        evolution.includeAmbiguousPresence = options.includeAmbiguousPresence;
        const AncestralOverlay overlay = AncestralInference::infer(graph, tree, evolution);
        const std::vector<AncestorChildAlignmentBlock> projections =
            AncestralInference::projectToIngroupChildren(graph, overlay);
        ancestorBlocks = overlay.blocks.size(); childBlocks = projections.size();
        writeAtomic(prefix + ".ancestor.calls.tsv",
                    [&](std::ostream &out) { AncestralWriters::writeCalls(overlay, out); });
        writeAtomic(prefix + ".ancestor.blocks.fa",
                    [&](std::ostream &out) { AncestralWriters::writeFasta(overlay, out); });
        writeAtomic(prefix + ".ancestor.children.maf", [&](std::ostream &out) {
            AncestralWriters::writeChildMaf(overlay, projections, out);
        });
        writeAtomic(prefix + ".ancestor.child-map.tsv", [&](std::ostream &out) {
            AncestralWriters::writeChildMap(projections, out);
        });
        if (!options.blockProjectionsPath.empty()) {
            AncestralAdjacencyConfig adjacencyConfig;
            adjacencyConfig.targetNode = options.ancestorNode;
            const AncestralAdjacencyReport adjacency =
                AncestralAdjacencyInference::infer(
                    graph, tree, structuralProjections, adjacencyConfig);
            adjacencyCalls = adjacency.calls.size();
            adjacencyIssues = adjacency.issues.size();
            for (const AncestralAdjacencyCall &call : adjacency.calls) {
                supportedAdjacencies += call.supportedWithoutLocalConflict;
            }
            writeAtomic(prefix + ".ancestor.block-projections.tsv",
                        [&](std::ostream &out) {
                AncestralAdjacencyWriters::writeProjections(
                    structuralProjections, out);
            });
            writeAtomic(prefix + ".ancestor.adjacencies.tsv",
                        [&](std::ostream &out) {
                AncestralAdjacencyWriters::writeCalls(adjacency, out);
            });
            writeAtomic(prefix + ".ancestor.adjacency-issues.tsv",
                        [&](std::ostream &out) {
                AncestralAdjacencyWriters::writeIssues(adjacency, out);
            });
        }
    }
    writeAtomic(prefix + ".qc.tsv", [&](std::ostream &out) {
        out << "metric\tvalue\n"
            << "core_hash\t" << graph.coreHash() << '\n'
            << "taxa\t" << taxa.records.size() << '\n'
            << "pairwise_runs\t" << pairs.records.size() << '\n'
            << "residues\t" << evidence.residues.size() << '\n'
            << "homology_edges\t" << evidence.homologies.size() << '\n'
            << "aligned_absence_observations\t" << evidence.alignedAbsences.size() << '\n'
            << "fasta_validation_enabled\t"
            << (options.validateFastaEvidence ? 1 : 0) << '\n'
            << "fasta_files_scanned\t" << fastaValidation.fastaFilesScanned << '\n'
            << "fasta_sources_validated\t" << fastaValidation.sourcesValidated << '\n'
            << "fasta_residues_validated\t"
            << fastaValidation.residueObservationsValidated << '\n'
            << "sites\t" << graph.sites().size() << '\n'
            << "conflicts\t" << graph.conflicts().size() << '\n'
            << "repair_candidates\t" << repairAudits.size() << '\n'
            << "repairs_applied\t";
        std::size_t repairsApplied = 0;
        for (const GraphRepairAudit &audit : repairAudits)
            repairsApplied += audit.disposition == "localized_poa_applied";
        out << repairsApplied << '\n'
            << "copy_status\t" << copyStatusName(copies.status) << '\n'
            << "extant_maf_blocks\t" << mafResult.blocksWritten << '\n'
            << "extant_maf_omissions\t" << mafResult.omissions.size() << '\n'
            << "ancestor_blocks\t" << ancestorBlocks << '\n'
            << "ancestor_child_blocks\t" << childBlocks << '\n'
            << "candidate_ancestral_adjacencies\t" << adjacencyCalls << '\n'
            << "supported_conflict_free_adjacencies\t" << supportedAdjacencies << '\n'
            << "adjacency_issues\t" << adjacencyIssues << '\n';
    });
    std::cout << "trioGraphAli completed; immutable core " << graph.coreHash()
              << ", sites " << graph.sites().size() << ", conflicts "
              << graph.conflicts().size() << '\n';
    return 0;
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave

int trioGraphAli(int argc, char **argv) {
    try {
        return anchorwave::trio::run(anchorwave::trio::parseOptions(argc, argv));
    } catch (const std::runtime_error &error) {
        if (std::string(error.what()) == "__help__") return 0;
        std::cerr << "trioGraphAli: " << error.what() << '\n';
    } catch (const std::exception &error) {
        std::cerr << "trioGraphAli: " << error.what() << '\n';
    }
    std::cerr << anchorwave::trio::usageText();
    return 1;
}
