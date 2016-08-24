UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
CC = g++
NORMAL_CFLAGS = -O2 -funroll-loops -w -g -fstack-protector-strong -Wformat -Werror=format-security
HDF5FLAG = -I/usr/include/hdf5/serial -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -Wdate-time -D_FORTIFY_SOURCE=2
HDF5LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a \
	/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a \
	-Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/seria
PNGLIBS = -lpng
OGL_LIBS = -lGL -lGLU -lglut -L/usr/X11R6/lib
endif


ifeq ($(UNAME), Darwin)
# do something Solaris-y
CC = g++
NORMAL_CFLAGS = -O3 -funroll-loops -w -g -fstack-protector-strong -Wformat -Werror=format-security
HDF5FLAG = -I/opt/local/include -I/opt/local/include -I/usr/include -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -Wdate-time -D_FORTIFY_SOURCE=2
HDF5LIBS = -L/opt/local/lib /opt/local/lib/libhdf5_hl.a /opt/local/lib/libhdf5.a -L/opt/local/lib -Wl,-headerpad_max_install_names -lpthread -lz -ldl -lm
PNGLIBS = -lpng
OGL_LIBS = -framework OpenGL -framework GLUT
endif

RM = rm -f

OUTPUT_EXTRACT_VIEW = lf_extract_view
OUTPUT_CORRES = lf_depth_v2
OUTPUT = $(OUTPUT_CORRES)

SRCS_COMMON = hdf5_lib/hdf5_tools.c png_lib/libPNG.c
OBJS_COMMON = $(subst .c,.o,$(notdir $(SRCS_COMMON)))

SRCS_CORRES = lf_depth_v2.c
OBJS_CORRES  = $(subst .c,.o,$(notdir $(SRCS_CORRES)))

SRCS_EXTRACT_VIEW = lf_extract_view.c
OBJS_EXTRACT_VIEW = $(subst .c,.o,$(notdir $(SRCS_EXTRACT_VIEW)))


#$(notdir $(FILE))
.PHONY: all clean

all: $(OUTPUT)

$(OUTPUT_CORRES): $(OBJS_COMMON) $(OBJS_CORRES)
	$(CC) $(NORMAL_CFLAGS) $(HDF5FLAG) -o $(OUTPUT_CORRES) $(OBJS_COMMON) $(OBJS_CORRES) $(HDF5LIBS) $(PNGLIBS) $(OGL_LIBS)

$(OBJS_COMMON): $(SRCS_COMMON)
	$(CC) $(NORMAL_CFLAGS) $(HDF5FLAG) -c $(SRCS_COMMON) $(HDF5LIBS) $(PNGLIBS)

$(OBJS_CORRES): $(SRCS_CORRES)
	$(CC) $(NORMAL_CFLAGS) $(HDF5FLAG) -c $(SRCS_CORRES) $(HDF5LIBS) $(PNGLIBS)


$(OUTPUT_EXTRACT_VIEW): $(OBJS_COMMON) $(OBJS_EXTRACT_VIEW)
	$(CC) $(NORMAL_CFLAGS) $(HDF5FLAG) -o $(OUTPUT_EXTRACT_VIEW) $(OBJS_EXTRACT_VIEW) $(OBJS_COMMON) $(HDF5LIBS) $(PNGLIBS)

$(OBJS_EXTRACT_VIEW): $(SRCS_EXTRACT_VIEW)
	$(CC) $(NORMAL_CFLAGS) $(HDF5FLAG) -c $(SRCS_EXTRACT_VIEW) $(HDF5LIBS) $(PNGLIBS)

clean:
	$(RM) $(OBJS) $(OBJS_EXTRACT_VIEW)
	$(RM) *.a
	$(RM) $(OUTPUT)
