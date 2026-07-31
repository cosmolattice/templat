# TempLat CI image -- incremental refresh.
#
# The cheap path, used by ci-image.yml on every push to main: start from the
# image already published for this config and just move its worktree to the new
# commit and rebuild what changed -- a few minutes instead of the ~11-14 min a
# rebuild from templat-ci-base.Dockerfile costs, and it skips the apt and
# FetchContent work entirely.
#
# Each refresh appends a Docker layer holding the objects and binaries that were
# rebuilt, so the image grows over time -- a merge touching a core header can add
# a few hundred MB. That growth is bounded by the weekly full rebuild in
# ci-image.yml, which is therefore load-bearing and not merely hygiene.
#
#   docker buildx build -t ghcr.io/cosmolattice/templat-ci:serial \
#       -f containers/ci/templat-ci-update.Dockerfile \
#       --build-arg base=ghcr.io/cosmolattice/templat-ci:serial \
#       --build-arg ref=<sha> --build-arg threads=8 .

ARG base=ghcr.io/cosmolattice/templat-ci:serial
FROM ${base}

ARG branch=main
ARG ref=main
ARG threads=4

# Fetch the branch first, then check out the commit: `git fetch origin <sha>`
# with a refspec is not portable, but once the branch is fetched the sha is
# reachable locally. `git checkout` rewrites only files whose content differs,
# which is exactly what keeps unchanged headers' mtimes -- and so make's
# incremental decisions -- intact. Re-running cmake picks up tests added or
# removed by the GLOB_RECURSE in tests/CMakeLists.txt; the cache already holds
# every -D flag, so none need repeating here.
RUN git -C /opt/templat/src fetch --no-tags --force origin "${branch}" \
    && git -C /opt/templat/src checkout --force "${ref}" \
    && git -C /opt/templat/src clean -xdff \
    && cmake -S /opt/templat/src -B /opt/templat/build \
    && ccache --zero-stats \
    && cmake --build /opt/templat/build -j "${threads}" \
    && ccache --show-stats

WORKDIR /opt/templat/build
