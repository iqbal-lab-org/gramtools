#include <gzstream.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <regex>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

// trim from start (in place)
static inline void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}


struct GenomicRead
{
	std::string *name;
	std::string *seq;
	std::string *qual;

	friend std::ostream& operator<< ( std::ostream & output, GenomicRead const & that )
	{ 
		return output << "[" << *(that.name) << "](" << *(that.seq) << ")";  
	}

	public:
	GenomicRead(){}

	GenomicRead (std::string n, std::string s,std::string q) 
	{
		name=new std::string(n);
		seq=new std::string(s);
		qual=new std::string(q);
	}

	~GenomicRead()
	{
		delete name;
		delete seq;
		delete qual;
	}

	std::vector<std::string> kmers (int k)
	{
		int sz=seq->length()-(k-1);
		std::vector<std::string> v=std::vector<std::string>();
		for (int i=0;i<sz;i++)
		{
			v.push_back(seq->substr(i,k));
		}
		return v;
	}

};


class SeqRead
{

	public:

		class EndOfFile : public std::runtime_error {
			public:
				EndOfFile( void ) : std::runtime_error("EndOfFile")
			{ }
		};

		class WrongInput : public std::runtime_error {
			public:
				WrongInput( void ) : std::runtime_error("WrongInput")
			{ }
		};

		class WrongFormat : public std::runtime_error {
			public:
				WrongFormat( void ) : std::runtime_error("WrongFormat")
			{ }
		};

		SeqRead (const char *fileinput)
		{
			filetype=0;
			gz=false;

			std::string file(fileinput);
			std::regex re_fastq(".*\\.(fastq|fq)(\\.gz)?$",std::regex::icase);
			std::regex re_fasta(".*\\.(fasta|fa)(\\.gz)?$",std::regex::icase);
			std::regex re_sam(".*\\.sam$",std::regex::icase);
			std::regex re_bam(".*\\.bam$",std::regex::icase);
			std::regex re_gz(".*\\.gz$",std::regex::icase);

			if (std::regex_match(file,re_sam))
				filetype=SAM;
			else if (std::regex_match(file,re_bam))
				filetype=BAM;
			else
			{
				if (std::regex_match(file,re_fastq))
					filetype=FASTQ;
				else if (std::regex_match(file,re_fasta))
					filetype=FASTA;
				else
					throw WrongInput();


				if (std::regex_match(file,re_gz))
					gz=true;
			}

			switch (filetype)
			{
				case SAM:
					break;
				case BAM:
					break;
				case FASTQ:
					nextproxy=&SeqRead::nextFQ;
					break;
				case FASTA:
					nextproxy=&SeqRead::nextFASTA;
					break;
			}

			if (filetype==FASTQ || filetype==FASTA)
			{
				if (gz)
				{
					inputptr=new igzstream();
					((igzstream*)inputptr)->open(fileinput);
				}
				else
				{
					inputptr=new std::ifstream();
					((std::ifstream*)inputptr)->open(fileinput);
				}
			}

		}

		~SeqRead()
		{
			if (filetype==FASTQ || filetype==FASTA)
			{ delete inputptr;}
		}


		class SeqIterator
		{

			public:
			SeqIterator (SeqRead* rd,long position=0) 
			{
				reader=rd;
				pos=position;
				if (position>=0)
				{
					try{ gr=reader->next();}
					catch (SeqRead::EndOfFile &e  )
					{pos=-1;}
				}
			}

			SeqIterator& operator++()
			{
				if (pos!=-1)
				{
					try{ gr=reader->next();}
					catch (SeqRead::EndOfFile &e  )
					{pos=-1;}
				}
				return *this;
			}

			bool operator==(const SeqIterator& rhs)
			{
				return rhs.pos==pos;
			}

			bool operator!=(const SeqIterator& rhs)
			{
				return rhs.pos!=pos;
			}

			GenomicRead* operator*()
			{
				return gr;
			}

			GenomicRead *gr;
			long pos;
			SeqRead *reader;
		};

		SeqIterator begin()
		{
			return SeqIterator(this,0);
		}

		SeqIterator end()
		{
			return SeqIterator(this,-1);
		}

		GenomicRead* next()
		{
			return (*this.*nextproxy)();
		}

		GenomicRead* nextFQ()
		{
			std::string n;
			std::string s;
			std::string q;

			if (!std::getline(*inputptr,n) || !std::getline(*inputptr,s) || !std::getline(*inputptr,q) || !std::getline(*inputptr,q))
				throw EndOfFile();


			return new GenomicRead(n,s,q);
		}

		GenomicRead* nextFASTA()
		{
			std::string s="";
			std::string n="";


			while (lastline.length()<1)
				if (!std::getline(*inputptr,lastline)) throw EndOfFile();

			if (lastline[0]!='>')
				throw WrongFormat();
			n=lastline.erase(0,1);

			while (1)
			{
				if (!std::getline(*inputptr,lastline))
					break;
				if (lastline.length() && lastline[0]!='>')
				{
					trim(lastline);
					s+=lastline;
				}
				else if (lastline.length() && lastline[0]=='>')
					break;
			}


			return new GenomicRead(n,s,"");
		}

		private:
		std::istream *inputptr;
		std::string lastline;

		GenomicRead* (SeqRead::*nextproxy)(void);

		bool gz;
		int filetype;
		static int const FASTA=1;
		static int const FASTQ=2;
		static int const BAM=3;
		static int const SAM=4;
};



int main (int argc, const char **argv)
{
	int a;
	SeqRead sq = SeqRead(argv[1]);
	for (auto rd: sq)
	{
//		std::cout << (*rd) << std::endl;
		for (auto kmer : rd->kmers(16))
		{
//			std::cout << kmer << std::endl;
		}
		delete rd;
	}
//	std::cout << "enter a\n";
//	std::cin >> a;
}
