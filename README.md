
```
bindgen PoaDraft.h \
  -o ../src/xpoa_sys.rs \
  --allowlist-function "PoaDraft.*" \
  -- \
  -x c++ \
  -std=c++17 \
  -I. \
  -I/usr/include/c++/13 \
  -I/usr/include/x86_64-linux-gnu/c++/13
```