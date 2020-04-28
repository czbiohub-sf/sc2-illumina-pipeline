#!/bin/bash
set -euxo pipefail
{
    function usage {
        echo "script usage: $(basename $0) [-r] [-o]

Arguments:
-r: Use nextflow -resume parameter
-o: Overwrite previous version-controlled test results
" >&2
    }

    RESUME=""
    OVERWRITE=""
    while getopts 'hro' opt; do
        case "$opt" in
            r) RESUME="-resume" ;;
            o) OVERWRITE=1 ;;
            h) usage
               exit 0;;
        esac
    done

    shift "$((OPTIND-1))"
    TEST_OUT="results/test"
    TEST_NEXTSTRAIN_OUT="results/test_nextstrain"

    TEST_OUT_STATS="$TEST_OUT/call_consensus-stats/combined.stats.tsv"
    TEST_NEXTSTRAIN_OUT_STATS="$TEST_NEXTSTRAIN_OUT/run_analysis-stats/combined.stats.tsv"

    TEST_STATS="data/test/test_results/test.combined.stats.tsv"
    TEST_NEXTSTRAIN_STATS="data/test/test_results/test_nextstrain.combined.stats.tsv"

    rm -rf "$TEST_OUT" "$TEST_NEXTSTRAIN_OUT"

    nextflow run call_consensus.nf "$RESUME" -profile test,docker
    nextflow run run_analysis.nf "$RESUME" -profile test_nextstrain,docker

    diff "$TEST_OUT_STATS" "$TEST_STATS" || true
    diff "$TEST_NEXTSTRAIN_OUT_STATS" "$TEST_NEXTSTRAIN_STATS" || true

    if [ -n "$OVERWRITE" ] ; then
        cp "$TEST_OUT_STATS" "$TEST_STATS"
        cp "$TEST_NEXTSTRAIN_OUT_STATS" "$TEST_NEXTSTRAIN_STATS"
    fi

    exit
}
