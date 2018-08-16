CC = mpic++
ARGS = -O3

all: cdice

cdice: CDICEMain.cpp CDICEInit.cpp CDICE.cpp CDICE_doeclim.cpp
	${CC} ${ARGS} -o a.exe CDICEMain.cpp CDICEInit.cpp CDICE.cpp CDICE_doeclim.cpp

clean:
	rm a.exe

