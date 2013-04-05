PROJECT(poisson)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# comment if you only going to write fortran code
ENABLE_LANGUAGE(C)

# needed no matter since the fst is written in fortran
ENABLE_LANGUAGE(Fortran)

ADD_EXECUTABLE(poisson poisson.c fst.f)
SET_TARGET_PROPERTIES(poisson PROPERTIES LINKER_LANGUAGE C)

ADD_EXECUTABLE(poisson_f poisson.f fst.f)
