#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <regex>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include "seq_file.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>



struct GenomicRead
{
	char *name;
	char *seq;
	char *qual;

	friend std::ostream& operator<< ( std::ostream & output, GenomicRead const & that )
	{ 
		return output << "[" << that.name << "](" << that.seq << ")";  
	}

	public:
	GenomicRead(){}

	std::vector<std::string> kmers (int k)
	{
		std::string tmp=std::string(seq);
		int sz=tmp.length()-(k-1);

		std::vector<std::string> v=std::vector<std::string>();
		for (int i=0;i<sz;i++)
		{
			v.push_back(tmp.substr(i,k));
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
			read = seq_read_new();
			file = seq_open(fileinput);
			gr=new GenomicRead();
		}

		~SeqRead()
		{
			seq_close(file);
			seq_read_free(read);
			free(gr);
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
			if (seq_read(file, read) > 0)
			{
				gr->name=read->name.b;
				gr->seq=read->seq.b;
				gr->qual=read->qual.b;
			}
			else
			{throw EndOfFile(); }
			return gr;
		}


		private:
		read_t *read;
		seq_file_t *file;
		GenomicRead *gr;
};



//int main (int argc, const char **argv)
//{
//	int a;
//	SeqRead sq = SeqRead(argv[1]);
//	for (auto rd: sq)
//	{
//		std::cout << (*rd) << std::endl;
//		for (auto kmer : rd->kmers(16))
//		{
//			std::cout << kmer << std::endl;
//		}
////		delete rd;
//	}
////	std::cout << "enter a\n";
////	std::cin >> a;
//}
