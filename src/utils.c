#include "methyl_extractor.h"
#include <stdarg.h>

void log_time(const char *format, ...) {
  time_t rawtime;
  struct tm *timeinfo;
  char time_str[64];
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", timeinfo);
  va_list args;
  va_start(args, format);
  fprintf(stderr, "[%s] ", time_str);
  vfprintf(stderr, format, args);
  va_end(args);
}

int make_directory(const char *path) {
  struct stat st = {0};
  if (stat(path, &st) == -1) {
    char *parent = strdup(path);
    char *slash = strrchr(parent, '/');
    if (slash && slash != parent) {
      *slash = '\0';
      make_directory(parent);
    }
    free(parent);
    return mkdir(path, 0700);
  }
  return 0;
}

uint8_t encode_nucleotide(char n) {
  switch (toupper(n)) {
  case 'A':
    return TNC_A;
  case 'C':
    return TNC_C;
  case 'G':
    return TNC_G;
  case 'T':
    return TNC_T;
  default:
    return TNC_N;
  }
}

char decode_nucleotide(uint8_t n) {
  switch (n) {
  case TNC_A:
    return 'A';
  case TNC_C:
    return 'C';
  case TNC_G:
    return 'G';
  case TNC_T:
    return 'T';
  default:
    return 'N';
  }
}

void decode_trinucleotide(uint8_t tnc, char *trinucl) {
  uint8_t n2 = (tnc / 5) % 4;
  uint8_t n3 = tnc % 5;
  trinucl[0] = 'C';
  trinucl[1] = decode_nucleotide(n2);
  trinucl[2] = decode_nucleotide(n3);
  trinucl[3] = '\0';
}

const char *get_context_string(int context) {
  if (context == CONTEXT_CPG)
    return "CG";
  else if (context == CONTEXT_CHG)
    return "CHG";
  else if (context == CONTEXT_CHH)
    return "CHH";
  return "???";
}

uint8_t encode_trinucleotide_context(const char *chr_seq, int pos, int chr_len,
                                     char strand) {
  int direction = (strand == '+') ? 1 : -1;
  uint8_t rv = 0;
  char base;

  if ((direction > 0 && pos + 2 >= chr_len) || (direction < 0 && pos <= 1))
    rv = 4;
  else {
    base = chr_seq[pos + 2 * direction];
    if (direction < 0)
      base = (base == 'A')   ? 'T'
             : (base == 'T') ? 'A'
             : (base == 'C') ? 'G'
             : (base == 'G') ? 'C'
                             : 'N';
    switch (toupper(base)) {
    case 'A':
      rv = 0;
      break;
    case 'C':
      rv = 1;
      break;
    case 'G':
      rv = 2;
      break;
    case 'T':
      rv = 3;
      break;
    default:
      rv = 4;
      break;
    }
  }

  if ((direction > 0 && pos + 1 >= chr_len) || (direction < 0 && pos == 0))
    rv += 20;
  else {
    base = chr_seq[pos + direction];
    if (direction < 0)
      base = (base == 'A')   ? 'T'
             : (base == 'T') ? 'A'
             : (base == 'C') ? 'G'
             : (base == 'G') ? 'C'
                             : 'N';
    switch (toupper(base)) {
    case 'A':
      rv += 0;
      break;
    case 'C':
      rv += 5;
      break;
    case 'G':
      rv += 10;
      break;
    case 'T':
      rv += 15;
      break;
    default:
      rv += 20;
      break;
    }
  }
  return rv;
}

int isCpG(const char *seq, int pos, int seqlen) {
  if (pos >= seqlen)
    return 0;
  if (toupper(seq[pos]) == 'C' && pos + 1 < seqlen &&
      toupper(seq[pos + 1]) == 'G')
    return 1;
  if (toupper(seq[pos]) == 'G' && pos > 0 && toupper(seq[pos - 1]) == 'C')
    return -1;
  return 0;
}

int isCHG(const char *seq, int pos, int seqlen) {
  if (pos + 2 >= seqlen)
    return 0;
  if (toupper(seq[pos]) == 'C' && toupper(seq[pos + 2]) == 'G')
    return 1;
  if (pos >= 2 && toupper(seq[pos]) == 'G' && toupper(seq[pos - 2]) == 'C')
    return -1;
  return 0;
}

int isCHH(const char *seq, int pos, int seqlen) {
  if (toupper(seq[pos]) == 'C')
    return 1;
  if (toupper(seq[pos]) == 'G')
    return -1;
  return 0;
}

int get_context(const char *chr_seq, int chr_len, int pos, int8_t *strand_ctx,
                uint8_t *tnc, int keep_chg, int keep_chh) {
  int context = 0;
  if (isCpG(chr_seq, pos, chr_len))
    context = CONTEXT_CPG;
  else if (keep_chg && isCHG(chr_seq, pos, chr_len))
    context = CONTEXT_CHG;
  else if (keep_chh && isCHH(chr_seq, pos, chr_len))
    context = CONTEXT_CHH;
  else
    return 0;

  char c2 = toupper(chr_seq[pos]);
  char strand = (c2 == 'C') ? '+' : '-';
  *tnc = encode_trinucleotide_context(chr_seq, pos, chr_len, strand);
  *strand_ctx = (strand == '+') ? context : -context;
  return context;
}

size_t count_methylation_sites(const char *chr_seq, uint32_t chr_len,
                               int keep_chg, int keep_chh) {
  size_t count = 0;
  for (uint32_t pos = 0; pos < chr_len; pos++) {
    if (isCpG(chr_seq, pos, chr_len) ||
        (keep_chg && isCHG(chr_seq, pos, chr_len)) ||
        (keep_chh && isCHH(chr_seq, pos, chr_len)))
      count++;
  }
  return count;
}

void initialize_buffer(MethylRecord *buffer, size_t site_count,
                       const char *chr_seq, uint32_t chr_len, int keep_chg,
                       int keep_chh) {
  size_t idx = 0;
  for (uint32_t pos = 0; pos < chr_len && idx < site_count; pos++) {
    int8_t strand_ctx;
    uint8_t tnc_val;
    int ctx = get_context(chr_seq, chr_len, pos, &strand_ctx, &tnc_val,
                          keep_chg, keep_chh);
    if (ctx) {
      buffer[idx].pos = pos + 1;
      buffer[idx].mC = buffer[idx].uC = 0;
      buffer[idx].tnc.tnc = tnc_val;
      buffer[idx].tnc.context = ctx;
      buffer[idx].tnc.strand = (strand_ctx > 0) ? 0 : 1;
      idx++;
    }
  }
}
