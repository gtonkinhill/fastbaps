/* The MIT License

 Copyright (c) 2008 Genome Research Ltd (GRL).

 Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/* Last Modified: 12APR2009 */

/* De-macro'd by the ACE-team 18MAY2012 wattup! */

/*Converted into template classes by CTS 11DEC2014*/

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <cctype>
#include <string>
#include <cstdlib>
#include <unistd.h>

#if HAVE_ZLIB
#include <zlib.h>
class FunctorZlib
{
public:
  int operator()(gzFile file, void * buffer, unsigned int len)
  {
    return gzread(file, buffer, len);
  }
};
#endif


#if HAVE_BZIP2
#include <bzlib.h>
class FunctorBZlib2
{
public:
  int operator()(BZFILE* file, void * buffer, int len)
  {
    return BZ2_bzread(file, buffer, len );
  }
};
#endif

class FunctorRead
{
public:
  ssize_t operator()(int fd, void *buf, size_t count)
  {
    return read(fd, buf, count);
  }
};

class kseq
{
public:
  kseq() : last_char(0) {};
  ~kseq() = default;;
  std::string name;
  std::string comment;
  std::string seq;
  std::string qual;
  int last_char;
};


template<class ret_t, class ReadFunction>
class kstream
{
public:
  kstream(ret_t f, ReadFunction rf) : bufferSize(4096), f(f), is_eof(0), begin(0), end(0), readfunc(rf) {
    buf = (char*) malloc(bufferSize);
  }

  ~kstream()
  {
    free(buf);
  }

  int read(kseq& seq)
  {
    int c;
    if (seq.last_char == 0)
    {
      while ((c = this->getc()) != -1 && c != '>' && c != '@');
      if (c == -1)
        return -1;
      seq.last_char = c;
    }
    seq.comment.clear();
    seq.seq.clear();
    seq.qual.clear();

    if (this->getuntil(0, seq.name, &c) < 0)
      return -1;
    if (c != '\n')
      this->getuntil( '\n', seq.comment, 0);
    while ((c = this->getc()) != -1 && c != '>' && c != '+' && c != '@')
    {
      if (isgraph(c))
      {
        seq.seq += (char)c;
      }
    }
    if (c == '>' || c == '@')
      seq.last_char = c;

    if (c != '+')
      return (int)seq.seq.length();


    while ((c = this->getc()) != -1 && c != '\n');

    if (c == -1)
      return -2;
    while ((c = this->getc()) != -1 && seq.qual.length() < seq.seq.length()) {
      if (c >= 33 && c <= 127)
        seq.qual += (char)c;
    }
    seq.last_char = 0;
    if (seq.seq.length() != seq.qual.length())
      return -2;
    return (int)seq.seq.length();
  }

private:
  int getc()
  {
    if (this->is_eof && this->begin >= this->end)
      return -1;
    if (this->begin >= this->end)
    {
      this->begin = 0;
      this->end = this->readfunc(this->f, this->buf, bufferSize);
      if (this->end < bufferSize)
        this->is_eof = 1;
      if (this->end == 0)
        return -1;
    }
    return (int)this->buf[this->begin++];
  }

  int getuntil(int delimiter, std::string &str, int *dret)
  {
    if (dret)
      *dret = 0;
    if (!str.empty()) {
      str.clear();
    }

    if (this->begin >= this->end && this->is_eof)
      return -1;
    for (;;)
    {
      int i;
      if (this->begin >= this->end)
      {
        if (!this->is_eof)
        {
          this->begin = 0;
          this->end = this->readfunc(this->f, this->buf, bufferSize);
          if (this->end < bufferSize)
            this->is_eof = 1;
          if (this->end == 0)
            break;
        }
        else
          break;
      }
      if (delimiter > 1)
      {
        for (i = this->begin; i < this->end; ++i)
        {
          if (this->buf[i] == delimiter)
            break;
        }
      }
      else if (delimiter == 0)
      {
        for (i = this->begin; i < this->end; ++i)
        {
          if (isspace(this->buf[i]))
            break;
        }
      }
      else if (delimiter == 1)
      {
        for (i = this->begin; i < this->end; ++i){
          if (isspace(this->buf[i]) && this->buf[i] != ' ')
            break;
        }
      }
      else i = 0;

      str.append(this->buf + this->begin, static_cast<unsigned long>(i - this->begin));
      this->begin = i + 1;
      if (i < this->end)
      {
        if (dret)
          *dret = this->buf[i];
        break;
      }
    }
    return (int)str.length();
  }

  char *buf;
  int begin;
  int end;
  int is_eof;
  ret_t f;
  ReadFunction readfunc;
  const unsigned int bufferSize;
};

#endif
