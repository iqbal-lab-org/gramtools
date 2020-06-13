#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <iostream>
#include <locale>
#include <regex>
#include <vector>
#include "seq_file.h"

#ifndef GRAMTOOLS_SEQREAD_HPP
#define GRAMTOOLS_SEQREAD_HPP

struct GenomicRead {
  std::string name;
  std::string seq;
  std::string qual;

  friend std::ostream &operator<<(std::ostream &output,
                                  GenomicRead const &that) {
    return output << "[" << that.name << "](" << that.seq << ")";
  }

 public:
  GenomicRead() {}
  GenomicRead(std::string const &name, std::string const &seq,
              std::string const &qual)
      : name(name), seq(seq), qual(qual) {}

  std::vector<std::string> kmers(int k) {
    std::string tmp = std::string(seq);
    int sz = tmp.length() - (k - 1);

    std::vector<std::string> v = std::vector<std::string>();
    for (int i = 0; i < sz; i++) {
      v.push_back(tmp.substr(i, k));
    }
    return v;
  }
};

using GenomicRead_vector = std::vector<GenomicRead>;

class AbstractGenomicReadIterator {
 public:
  virtual AbstractGenomicReadIterator &operator++() = 0;

  /**
   * Use this as end condition in for/while loops
   */
  bool has_more_reads() { return pos != -1; }

  GenomicRead *operator*() { return gr; }

 protected:
  long pos;
  GenomicRead *gr;
};

/**
 * Takes existing specified reads and makes an iterator consistent with
 * SeqIterator below, which works on files.
 */
class GenomicReadIterator : public AbstractGenomicReadIterator {
 public:
  GenomicReadIterator(GenomicRead_vector const &input_reads)
      : reads(input_reads) {
    num_reads = reads.size();
    assert(num_reads > 0);
    pos = 0;
    gr = &reads[pos];
  }

  AbstractGenomicReadIterator &operator++() override {
    if (pos != -1) {
      if (pos < num_reads - 1) {
        pos++;
        gr = &reads[pos];
      } else if (pos == num_reads - 1)
        pos = -1;
    }
    return *this;
  }

 private:
  GenomicRead_vector reads;
  std::size_t num_reads;
};

class SeqRead {
 public:
  class EndOfFile : public std::runtime_error {
   public:
    EndOfFile(void) : std::runtime_error("EndOfFile") {}
  };

  class WrongInput : public std::runtime_error {
   public:
    WrongInput(void) : std::runtime_error("WrongInput") {}
  };

  class WrongFormat : public std::runtime_error {
   public:
    WrongFormat(void) : std::runtime_error("WrongFormat") {}
  };

  SeqRead(const char *fileinput) {
    read = seq_read_new();
    file = seq_open(fileinput);
    if (file == NULL) {
      seq_read_free(read);
      printf("Unable to open %s\n", fileinput);
      exit(1);
    } else {
      gr = new GenomicRead();
    }
  }

  ~SeqRead() {
    seq_close(file);
    seq_read_free(read);
    delete gr;
  }

  class SeqIterator : public AbstractGenomicReadIterator {
   public:
    SeqIterator(SeqRead *rd, long position = 0) {
      reader = rd;
      pos = position;
      if (position >= 0) {
        try {
          gr = reader->next();
        } catch (SeqRead::EndOfFile &e) {
          pos = -1;
        }
      }
    }

    SeqIterator &operator++() override {
      if (pos != -1) {
        try {
          gr = reader->next();
        } catch (SeqRead::EndOfFile &e) {
          pos = -1;
        }
      }
      return *this;
    }

    bool operator==(const SeqIterator &rhs) { return rhs.pos == pos; }

    bool operator!=(const SeqIterator &rhs) { return rhs.pos != pos; }

    SeqRead *reader;
  };

  SeqIterator begin() { return SeqIterator(this, 0); }

  SeqIterator end() { return SeqIterator(this, -1); }

  GenomicRead *next() {
    if (seq_read(file, read) > 0) {
      gr->name = read->name.b;
      gr->seq = read->seq.b;
      gr->qual = read->qual.b;
    } else {
      throw EndOfFile();
    }
    return gr;
  }

 private:
  read_t *read;
  seq_file_t *file;
  GenomicRead *gr;
};

#endif  // GRAMTOOLS_SEQREAD_HPP
