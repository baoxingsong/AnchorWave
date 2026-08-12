if(NOT DEFINED ANCHORWAVE_EXE OR NOT DEFINED FIXTURE_DIR OR NOT DEFINED OUTPUT_PREFIX)
  message(FATAL_ERROR "CLI smoke test variables are missing")
endif()

# The experimental command remains directly callable for regression testing,
# but it is deliberately absent from the public top-level help.
execute_process(
  COMMAND "${ANCHORWAVE_EXE}" --help
  RESULT_VARIABLE help_status
  OUTPUT_VARIABLE help_stdout
  ERROR_VARIABLE help_stderr)
if(NOT help_status EQUAL 0)
  message(FATAL_ERROR "anchorwave --help failed (${help_status})\n${help_stdout}\n${help_stderr}")
endif()
if(help_stdout MATCHES "trioGraphAli" OR help_stderr MATCHES "trioGraphAli")
  message(FATAL_ERROR "experimental trioGraphAli leaked into public help")
endif()

file(GLOB old_outputs "${OUTPUT_PREFIX}.*")
if(old_outputs)
  file(REMOVE ${old_outputs})
endif()

execute_process(
  COMMAND "${ANCHORWAVE_EXE}" trioGraphAli
          --taxa "${FIXTURE_DIR}/taxa.tsv"
          --pairwise-manifest "${FIXTURE_DIR}/pairwise.tsv"
          --species-tree "${FIXTURE_DIR}/species.nwk"
          --copy-relations "${FIXTURE_DIR}/copy_relations.tsv"
          --block-projections "${FIXTURE_DIR}/block_projections.tsv"
          --copy-mode strict
          --ancestor-node A
          --output-prefix "${OUTPUT_PREFIX}"
  RESULT_VARIABLE command_status
  OUTPUT_VARIABLE command_stdout
  ERROR_VARIABLE command_stderr)
if(NOT command_status EQUAL 0)
  message(FATAL_ERROR "trioGraphAli failed (${command_status})\n${command_stdout}\n${command_stderr}")
endif()

set(required_suffixes
    graph.gfa graph.meta.tsv extant.maf extant.maf.omissions.tsv
    copy-resolution.tsv discordance.tsv
    repair-audit.tsv qc.tsv ancestor.calls.tsv ancestor.blocks.fa
    ancestor.children.maf ancestor.child-map.tsv ancestor.block-projections.tsv
    ancestor.adjacencies.tsv ancestor.adjacency-issues.tsv)
foreach(suffix IN LISTS required_suffixes)
  set(path "${OUTPUT_PREFIX}.${suffix}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "expected output is missing: ${path}")
  endif()
  file(SIZE "${path}" output_size)
  if(output_size EQUAL 0)
    message(FATAL_ERROR "expected output is empty: ${path}")
  endif()
endforeach()

file(READ "${OUTPUT_PREFIX}.qc.tsv" qc)
if(NOT qc MATCHES "core_hash.core_")
  message(FATAL_ERROR "QC output lacks immutable core hash")
endif()
file(READ "${OUTPUT_PREFIX}.ancestor.children.maf" child_maf)
if(NOT child_maf MATCHES "s.A\\.ancestor_block_")
  message(FATAL_ERROR "ancestor-child MAF lacks ancestor row")
endif()
