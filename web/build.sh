#!/usr/bin/env bash
# Builds the site into dist/.   ./build.sh        -> web page build
#                                         ./build.sh test   -> node build + smoke test
#                                         ./build.sh serve  -> build, then serve dist/ on :8000 (or $PORT; next free if taken)
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v emcmake >/dev/null; then
  # shellcheck disable=SC1091
  source "${EMSDK:-$HOME/emsdk}/emsdk_env.sh" >/dev/null 2>&1
fi

case "${1:-}" in
  test)
    emcmake cmake -S . -B build-node -DHLSIM_ENVIRONMENT=node -DHLSIM_OUT_DIR="$PWD/build-node/out" >/dev/null
    cmake --build build-node -j
    node test/smoke.mjs "${@:2}"
    ;;
  *)
    emcmake cmake -S . -B build >/dev/null
    cmake --build build -j
    if [[ "${1:-}" == serve ]]; then
      port="${PORT:-8000}"
      # Skip ports already in use (e.g. a server left over from an earlier run).
      while ! python3 -c 'import socket,sys; s=socket.socket(); s.bind(("", int(sys.argv[1])))' "$port" 2>/dev/null; do
        echo "port $port in use, trying $((port + 1))" >&2
        port=$((port + 1))
      done
      echo "http://localhost:$port"
      exec python3 -m http.server -d dist "$port"
    fi
    ;;
esac
