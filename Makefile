CC = g++
ARGS = -O3

all: cdice

cdice: CDICEMain.cpp CDICEInit.cpp CDICE.cpp
	${CC} ${ARGS} -o cdice CDICEMain.cpp CDICEInit.cpp CDICE.cpp

clean:
	rm cdice

