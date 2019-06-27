noinst_LTLIBRARIES = libtest.la libfv_mapz.la
libfv_mapz_la_SOURCES = fv_mapz.F90
libfv_mapz_la_LIDADD = ./libtest.la

libtest_la_SOURCES = test.F90
include_HEADERS = fv_mapz_mod.mod test.mod

bin_PROGRAMS = fv_mapz
fv_mapz_SOURCES = main.F90
fv_mapz_LDADD = ./libfv_mapz.la ./libtest.la
