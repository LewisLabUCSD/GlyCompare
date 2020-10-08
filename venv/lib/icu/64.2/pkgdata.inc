GENCCODE_ASSEMBLY_TYPE=-a gcc-darwin
SO=dylib
SOBJ=dylib
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.64.2.dylib
COMPILE=/usr/local/miniconda/conda-bld/icu_1566795134219/_build_env/bin/clang -DU_ATTRIBUTE_DEPRECATED= -D_FORTIFY_SOURCE=2 -mmacosx-version-min=10.9 -I/Users/zhangyujie/anaconda3/envs/glycompare_env_2   -DU_HAVE_STRTOD_L=1 -DU_HAVE_XLOCALE_H=1  -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fPIE -fstack-protector-strong -O2 -pipe -I/Users/zhangyujie/anaconda3/envs/glycompare_env_2/include -fdebug-prefix-map=/usr/local/miniconda/conda-bld/icu_1566795134219/work=/usr/local/src/conda/icu-64.2 -fdebug-prefix-map=/Users/zhangyujie/anaconda3/envs/glycompare_env_2=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -fno-common -c
LIBFLAGS=-I/Users/zhangyujie/anaconda3/envs/glycompare_env_2/include  -dynamic
GENLIB=/usr/local/miniconda/conda-bld/icu_1566795134219/_build_env/bin/clang -dynamiclib -dynamic -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fPIE -fstack-protector-strong -O2 -pipe -I/Users/zhangyujie/anaconda3/envs/glycompare_env_2/include -fdebug-prefix-map=/usr/local/miniconda/conda-bld/icu_1566795134219/work=/usr/local/src/conda/icu-64.2 -fdebug-prefix-map=/Users/zhangyujie/anaconda3/envs/glycompare_env_2=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-pie -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/zhangyujie/anaconda3/envs/glycompare_env_2/lib -L/Users/zhangyujie/anaconda3/envs/glycompare_env_2/lib -Wl,-rpath,/Users/zhangyujie/anaconda3/envs/glycompare_env_2/lib 
LDICUDTFLAGS=
LD_SONAME=-Wl,-compatibility_version -Wl,64 -Wl,-current_version -Wl,64.2 -install_name /Users/zhangyujie/anaconda3/envs/glycompare_env_2/lib/
RPATH_FLAGS=
BIR_LDFLAGS=
AR=/usr/local/miniconda/conda-bld/icu_1566795134219/_build_env/bin/x86_64-apple-darwin13.4.0-ar
ARFLAGS=r -c
RANLIB=/usr/local/miniconda/conda-bld/icu_1566795134219/_build_env/bin/x86_64-apple-darwin13.4.0-ranlib
INSTALL_CMD=/usr/bin/install -c
