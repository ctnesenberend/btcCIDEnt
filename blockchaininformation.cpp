#include <iostream>
#include <fstream> 
#include <map>
#include <cstdint>
#include <string>
#include <cstring>
#include <vector> 
#include <zlib.h>
#include <time.h>
#include <cmath>

class datastream {
	// Basic class to consider the datafiles, with on each row integers (possibly except one)
	private:
	std::ifstream stream;
	std::string sdata;
	std::vector<uint64_t> data;
	int length;
	int skip;
	bool open;
	
	public:
	datastream(std::string filename, int l, int k = -1){
		openstream(filename,l,k);
	}
	
	datastream(char* filename, int l, int k = -1){
		openstream(filename,l,k);
	}
	
	datastream(){
		skip = -1;
		open = 0;
		length = 0;
	}
	
	void openstream(char* filename, int l, int k = -1){
		// opens file filename which has l columns (possibly skipping column k; the skip is to skip adresses which aren't integers)
		skip = k;
		length = l;
		open = 1;
		data.resize(l);
		stream.open(filename);
		if(!stream.is_open()){
			std::cout << "Failed to load file "<< filename << ", exiting now" << std::endl;
			exit(0);
		}		
	}
	
	void openstream(std::string filename, int l, int k = -1){
		// opens file filename which has l columns (possibly skipping column k; the skip is to skip adresses which aren't integers)		
		skip = k;
		length = l;
		open = 1;
		data.resize(l);
		stream.open(filename);
		if(!stream.is_open()){
			std::cout << "Failed to load file "<< filename << ", exiting now" << std::endl;
			exit(0);
		}		
	}
	void newl(){
		// fetches a new line of data from the file and stores it in the class 
		if(!open)
			return;
		if(!stream.good()){
			std::cout << "Trying to stream out of bounds, ending program" << std::endl;
			exit(1);
		}
		for (int i = 0; i< length; i++){
			if(moredata())
				stream >> sdata;
			if(i != skip)
				data[i] = std::stoull(sdata);
		}
	}
	uint64_t &operator[](int i){
		// fetch the ith element of the currently stored line
		if (!open || i >= length){
			std::cout << "Trying to stream out of bounds, ending program" << std::endl;
			exit(1);
		}
		return data[i];
	}
	
	bool moredata(){
		// checks if there's more data left
		return open && stream.peek()!=EOF;
		
	}
	
	bool is_open(){
		return open && stream.is_open();
	}
	
};

class infoStream{
	// class that evaluates the blockchain and calculates CID 
	private:

	std::map<uint64_t,uint64_t> money;
	std::vector<unsigned char> moneys;
	std::map<uint64_t,uint64_t> moneyvals;
	std::vector<unsigned char> compr;
	bool getentropy;
	bool getCID;
	uint64_t maxmoney;
	
	public:
	datastream in;
	datastream out;
	datastream block;
	datastream trans;
	uint64_t time;
	uint64_t folx;
	uint64_t nusers;
	long double de;
	long double entropyval;

	infoStream(char* instring, char* outstring, char* blockstring, char* transstring){
		// Opens stream based on inputs 
		

		
		// open streams as required 
		in.openstream(instring,6);
		out.openstream(outstring,4);
		block.openstream(blockstring,4,1); // the second column is an adress and therefor not an integer 
		trans.openstream(transstring,4);
		
		// prefetch the first line 
		in.newl();
		out.newl();
		block.newl();
		trans.newl();
		time = 0;
		
		// stuff for CID calculation
		nusers = 0;
		moneys.resize(8*20000000); // this is the array in which the to be compressed data is stored 
		compr.resize(compressBound(8*20000000)); // and this in which  the compressed data is stored 

		// stuff for real-time entropy calculation
		folx = 0;
		time = 0;
		de = 0; // this is -k * log(k)

	}
	

	void dotrans(uint64_t timestamp){
		// Calculates transactions up to and including the first block >= timestamp
		while(time <= timestamp){
			oneblock();
		}
	}
	
	bool oneblock(){
		// Applies transactions from exactly one block and returns whether there is more data 
		time = block[2];
		// apply the transactions from the block
		for (; trans[1] == block[0]; trans.newl()){
			// go do the + transactions of the transaction
			for(; trans[0] == out[0]; out.newl()){
				uint64_t money_cur = money[out[2]];
				uint64_t money_fut = money_cur + out[3];
				
				if(money_cur == 0)
					folx++;
				else {
					moneyvals[money_cur]--;
					if(moneyvals[money_cur] != 0)
						de += ((long double) moneyvals[money_cur])*logl(1+1/((long double) moneyvals[money_cur])) + logl((long double) moneyvals[money_cur]+1); // addition as de is always negative
				}
				
				if(moneyvals[money_cur] == 0)
					moneyvals.erase(money_cur);
				
				
				if(moneyvals[money_fut] != 0)  {
					de += -((long double)moneyvals[money_fut])*logl(1+1/((long double) moneyvals[money_fut])) - logl((long double) (moneyvals[money_fut]+1));
				}
				moneyvals[money_fut]++;
				
				money[out[2]] += out[3];

			}
			// do all the - transactions of the transaction
			for(; trans[0] == in[0]; in.newl()){
				uint64_t money_cur = money[in[4]];
				uint64_t money_fut = money_cur - in[5];
				
				moneyvals[money_cur]--;
				if(moneyvals[money_cur] != 0)
					de += ((long double)moneyvals[money_cur])*logl(1+1/((long double) moneyvals[money_cur])) + logl((long double) moneyvals[money_cur]+1);
				
				
				if(moneyvals[money_cur] == 0)
					moneyvals.erase(money_cur);
				
				if(money_fut == 0)
					folx--;
				else{
					if(moneyvals[money_fut] != 0)  
						de += -((long double)moneyvals[money_fut])*logl(1+1/((long double) moneyvals[money_fut])) - logl((long double) moneyvals[money_fut]+1);
					moneyvals[money_fut]++;
				}
				
				
				
				money[in[4]] -= in[5];
				// delete any elements that end up zero, important for memory concerns (no need to do ths with the + transactions obviously)
				if(money[in[4]]== 0) 
					money.erase(in[4]);
			}
		}
		if(block.moredata()){
			block.newl(); // prefetches the next block
			return 1;
		}
		else return 0;
	}

	void fetcharray(){
		// turns the map into a char array for compression purposes
		nusers = 0;
		for( const auto &m : money){
			memcpy(&moneys[8*nusers],&m.second,8);
			nusers++;
			if(nusers*8 >= moneys.size()){ // very basic dynamic array; no need for it to be more complex
				moneys.resize(moneys.size() + 8*20000000);
				compr.resize(compressBound(moneys.size()));
			}
		}
	}
		
	
	double CID(){
		// calculates the CID at the current timestamp 
		fetcharray();
		uLongf lengths = nusers*8; // *8 because a uint64_t is 8 bytes
		uLongf lengthc = compressBound(lengths);
		compress(&compr[0], &lengthc, (const Bytef*) &moneys[0], lengths);
		return ((double)lengthc)/((double) lengths);

	
	}
	double entropy(){
		// calculates the entropy at the current timestamp
		entropyval =  logl( (long double) folx) + de/folx; 
		return entropyval;
	}
	

};

int main(int argc, char** argv){
	std::string str_txin = std::string("txin.dat");
	std::string str_txout = std::string("txout.dat");
	std::string str_bh = std::string("bh.dat");
	std::string str_tx = std::string("tx.dat");
	
	if(argc >= 5){
		str_txin = std::string(argv[1]);
		str_txout = std::string(argv[2]);
		str_bh = std::string(argv[3]);
		str_tx = std::string(argv[4]);
	}
	
	infoStream strm(argv[1],argv[2],argv[3],argv[4]); // initialize streams 
	
	uint64_t begintime = 0;
	if(argc > 5)
		begintime = strtoull(argv[5],NULL,10); // Allows setting a starting time of calculating the CID after the start, to be able to do the full calculation over multiple runs
	if(argc == 2)
		begintime = strtoull(argv[1],NULL,10);
	
	std::string filename = std::string("blockchaininfo.dat");
	if(begintime != 0)
		filename = std::string("blockchaininfo")+std::to_string(begintime)+std::string(".dat");
	
	std::ofstream ofile(filename, std::ofstream::out | std::ofstream::trunc);
	
	// Evaluates all blocks
	ofile.precision(15);
	ofile << std::fixed; // have consistent notation of the floats 
	
	for(strm.oneblock(); strm.time < 1518086764; strm.oneblock()){ // loop over the blocks 
		std::cout << "\rEvaluating current block at: " << strm.time;
		if (strm.time > begintime){ // only evaluate blocks after the defined begintime 
			ofile << strm.time << //"\t" << strm.CID() << 
				"\t" << strm.entropy() << "\t" << strm.nusers << std::endl;
		}
	}
	std::cout << "\rDone!                                   " << std::endl;
	return 0;
}