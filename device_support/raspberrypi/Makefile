all: libAsyncStoreManager libBCM2835Interface
	echo Done.


libAsyncStoreManager: AsyncStoreManager.cpp
	gcc -g -O3 -fpic -shared -pthread -o libAsyncStoreManager.so AsyncStoreManager.cpp -I${MDSPLUS_DIR}/include -L${MDSPLUS_DIR}/lib -lm -lpthread -lMdsObjectsCppShr
	cp libAsyncStoreManager.so /usr/local/mdsplus/lib

libBCM2835Interface: BCM2835Interface.cpp
	gcc    -g -O3 -fpic -shared -pthread -o libBCM2835Interface.so BCM2835Interface.cpp -L${MDSPLUS_DIR}/lib -lbcm2835 -lrt -lAsyncStoreManager -I${MDSPLUS_DIR}/include
	#gcc    -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -fpic -fsigned-char -shared-libgcc -fno-strict-aliasing -pthread -shared -o libBCM2835Interface.so BCM2835Interface.cpp -lbcm2835
	sudo chown root:root /usr/local/mdsplus/bin/mdsdcl
	sudo chmod +s /usr/local/mdsplus/bin/mdsdcl
	cp libBCM2835Interface.so /usr/local/mdsplus/lib

