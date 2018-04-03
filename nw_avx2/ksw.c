//Bampetas Nikolaos 4511859
//Research program in interest of computer engineering lab TU Delft


#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
//#include <emmintrin.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>
#include "kseq.h"
#include "ksw.h"
#include "utils.h"
#include <immintrin.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif
//#include <zmmintrin.h>  avx512


#define MATCH 1
#define MISMATCH -1
#define GAPO  6
#define GAPE  1
#define SIZE 1000
#define NEG_INF (INT16_MIN/(int16_t)(2))
#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

typedef union __m256i_16 {
    __m256i m;
    int16_t v[16];
} __m256i_16_t;

typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp ()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}
static inline __m256i _mm256_insert_epi16_rpl(__m256i a, int16_t i, int imm) {
    __m256i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
    A.m = a;
    return A.v[imm];
}
struct _kswq_t {
   int qlen, slen;
   uint8_t shift, mdiff, max, size;
   __m256i *qp, *H0, *H1, *E, *Hmax;
};

const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 };
KSEQ_INIT(gzFile, err_gzread)
int qlen;
kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry);
kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int xtra,kswq_t **qry);
static inline void revseq(int l, uint8_t *s) ;
kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat);
kswr_t ksw_u8(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra);
kswr_t ksw_i16(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra);
 FILE *f;
int main(int argc,char* argv[]){
timestamp_t t0, t1,secs=0;

	int8_t mat[25];	
	int gapo = 6, gape = 1, xtra = KSW_XSTART;
       int minsc = 0;
	kseq_t *kst, *ksq;
	gzFile fpt, fpq;
	int i,k,j, sa=1,sb=4;
		
    f = fopen("test.txt", "w");
	if(argc!=3){
		printf("wrong number of arguments\n");
	}else{
	
		// initialize scoring matrix
   		for (i = k = 0; i < 4; ++i) {
     			for (j = 0; j < 4; ++j)
      				mat[k++] = i == j? sa : -sb;
      			mat[k++] = 0; // ambiguous base
		}
		for (j = 0; j < 5; ++j) 
			mat[k++] = 0;
		//-------------------------------------------
		//-----READ DNA DATA FROM FASTA FILES------------
		fpt = xzopen(argv[1], "r"); kst = kseq_init(fpt);
  		fpq = xzopen(argv[2], "r"); ksq = kseq_init(fpq);
		
		while (kseq_read(ksq) > 0 && kseq_read(kst) > 0) {
			kswr_t r;
			r.score = -1, r.te = -1, r.qe = -1, r.qb = -1, r.tb = -1;
			for (i = 0; i < (int)ksq->seq.l; ++i) ksq->seq.s[i] = seq_nt4_table[(int)ksq->seq.s[i]];
			for (i = 0; i < (int)kst->seq.l; ++i) kst->seq.s[i] = seq_nt4_table[(int)kst->seq.s[i]];
			qlen=ksq->seq.l;
			t0 = get_timestamp();
			r = ksw_align(ksq->seq.l, (uint8_t*)ksq->seq.s, kst->seq.l, (uint8_t*)kst->seq.s, 5, mat, gapo, gape, xtra, 0);
			t1 = get_timestamp();
			secs += (t1 - t0);
			//printf(" %lld usecs real time \n", secs);
			//if (r.score == -1 || r.te == -1 || r.qe == -1 || r.qb == -1 || r.tb == -1) {
        	//	fprintf(stderr, "AVX2 implementation not working\n");
        			//exit(EXIT_FAILURE);
      			//}
      			if (r.score >= minsc)
   			 printf("ref name=%s\tread name=%s\tscore=%d\n", kst->name.s, ksq->name.s, r.score);
		}
		
  	 kseq_destroy(kst); err_gzclose(fpt);
   	 kseq_destroy(ksq); err_gzclose(fpq);
fclose(f);
	}
	return 0;
}




kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry) {
   return ksw_align2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, xtra, qry);
}

kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int xtra,kswq_t **qry) {


   kswq_t *q;
   kswr_t r;
   kswr_t (*func)(kswq_t*, int, const uint8_t*, int, int, int, int, int);
   
   q = (qry && *qry) ? *qry : ksw_qinit((xtra & KSW_XBYTE) ? 1 : 2, qlen, query, m, mat);
   if (qry && *qry == 0)
      *qry = q;    
   func = q->size == 2 ? ksw_i16 : ksw_u8;   	
   r = func(q, tlen, target, o_del, e_del, o_ins, e_ins, xtra);
    if (qry == 0)
      free(q);
   if ((xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff)))
      return r;
   return r;
}

static inline void revseq(int l, uint8_t *s) {
   int i, t;
   for (i = 0; i < l >> 1; ++i)
      t = s[i], s[i] = s[l - 1 - i], s[l - 1 - i] = t;
	 
}

kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat) {
   kswq_t *q;
   int slen, a, tmp, p;
   size = size > 1 ? 2 : 1;
   p = 8 * (4 - size); // # values per __m256i (16 16bit or 32 8bit)
   slen = (qlen + p - 1) / p; // segmented length
   q = (kswq_t*) malloc(sizeof(kswq_t) + 256 + 32 * slen * (m + 4)); // a single block of memory
   q->qp = (__m256i*) (((size_t) q + sizeof(kswq_t) + 31) >> 5 << 5); // align memory   32bit boundary
   q->H0 = q->qp + slen * m;
   q->H1 = q->H0 + slen;
   q->E = q->H1 + slen;
   q->Hmax = q->E + slen;
   q->slen = slen;
   q->qlen = qlen;
   q->size = size;
   // compute shift
   tmp = m * m;
   for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
      if (mat[a] < (int8_t) q->shift)
         q->shift = mat[a];
      if (mat[a] > (int8_t) q->mdiff)
         q->mdiff = mat[a];
   }
   q->max = q->mdiff;
   q->shift = 256 - q->shift; // NB: q->shift is uint8_t
   q->mdiff += q->shift; // this is the difference between the min and max scores
   // An example: p=8, qlen=19, slen=3 and segmentation:
   //  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
   if (size == 1) {
      int8_t *t = (int8_t*) q->qp;
      for (a = 0; a < m; ++a) {
         int i, k, nlen = slen * p;
         const int8_t *ma = mat + a * m;
         for (i = 0; i < slen; ++i)
            for (k = i; k < nlen; k += slen) // p iterations
               *t++ = (k >= qlen ? 0 : ma[query[k]]) + q->shift;
      }
   } else {
      int16_t *t = (int16_t*) q->qp;
      for (a = 0; a < m; ++a) {
         int i, k, nlen = slen * p;
         const int8_t *ma = mat + a * m;
         for (i = 0; i < slen; ++i)
            for (k = i; k < nlen; k += slen) // p iterations
               *t++ = (k >= qlen ? NEG_INF : ma[query[k]]);
      }
   }
 
   return q;
}

kswr_t ksw_u8(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
   int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
   uint64_t *b;
   __m256i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
   kswr_t r;


#define __max_16(ret, xx) do { \
		(xx) = _mm256_max_epu8((xx), _mm256_srli_si256((xx), 16)); \
		(xx) = _mm256_max_epu8((xx), _mm256_srli_si256((xx), 8)); \
		(xx) = _mm256_max_epu8((xx), _mm256_srli_si256((xx), 4)); \
		(xx) = _mm256_max_epu8((xx), _mm256_srli_si256((xx), 2)); \
		(xx) = _mm256_max_epu8((xx), _mm256_srli_si256((xx), 1)); \
    	(ret) = _mm256_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

   // initialization
   r = g_defr;
   minsc = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
   endsc = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
   m_b = n_b = 0;
   b = 0;
   zero = _mm256_set1_epi32(0);
   oe_del = _mm256_set1_epi8(_o_del + _e_del);
   e_del = _mm256_set1_epi8(_e_del);
   oe_ins = _mm256_set1_epi8(_o_ins + _e_ins);
   e_ins = _mm256_set1_epi8(_e_ins);
   shift = _mm256_set1_epi8(q->shift);
   H0 = q->H0;
   H1 = q->H1;
   E = q->E;
   Hmax = q->Hmax;
   slen = q->slen;

   for (i = 0; i < slen; ++i) {
      _mm256_store_si256(E + i, zero);
      _mm256_store_si256(H0 + i, zero);
      _mm256_store_si256(Hmax + i, zero);
   }
   // the core loop
   for (i = 0; i < tlen; ++i) {
      int j, k, cmp, imax;
      __m256i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
      h = _mm256_load_si256(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
      h = _mm256_slli_si256(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
      for (j = 0; LIKELY(j < slen); ++j) {
         /* SW cells are computed in the following order:
          *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
          *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
          *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
          */

         // compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
         h = _mm256_adds_epu8(h, _mm256_load_si256(S + j));
         h = _mm256_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
         e = _mm256_load_si256(E + j); // e=E'(i,j)
         h = _mm256_max_epu8(h, e);
         h = _mm256_max_epu8(h, f); // h=H'(i,j)
         max = _mm256_max_epu8(max, h); // set max
         _mm256_store_si256(H1 + j, h); // save to H'(i,j)
         // now compute E'(i+1,j)
         e = _mm256_subs_epu8(e, e_del); // e=E'(i,j) - e_del
         t = _mm256_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
         e = _mm256_max_epu8(e, t); // e=E'(i+1,j)
         _mm256_store_si256(E + j, e); // save to E'(i+1,j)
         // now compute F'(i,j+1)
         f = _mm256_subs_epu8(f, e_ins);
         t = _mm256_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
         f = _mm256_max_epu8(f, t);
         // get H'(i-1,j) and prepare for the next j
         h = _mm256_load_si256(H0 + j); // h=H'(i-1,j)
      }
      // NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
      for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
         f = _mm256_slli_si256(f, 1);
         for (j = 0; LIKELY(j < slen); ++j) {
            h = _mm256_load_si256(H1 + j);
            h = _mm256_max_epu8(h, f); // h=H'(i,j)
            _mm256_store_si256(H1 + j, h);
            h = _mm256_subs_epu8(h, oe_ins);
            f = _mm256_subs_epu8(f, e_ins);
            cmp = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_subs_epu8(f, h), zero));
            if (UNLIKELY(cmp == 0xfffff))
               goto end_loop16;
         }
      }
      end_loop16:
      //int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
      __max_16(imax, max); // imax is the maximum number in max
      if (imax >= minsc) { // write the b array; this condition adds branching unfornately
         if (n_b == 0 || (int32_t) b[n_b - 1] + 1 != i) { // then append
            if (n_b == m_b) {
               m_b = m_b ? m_b << 1 : 8;
               b = (uint64_t*) realloc(b, 8 * m_b);
            }
            b[n_b++] = (uint64_t) imax << 32 | i;
         } else if ((int) (b[n_b - 1] >> 32) < imax)
            b[n_b - 1] = (uint64_t) imax << 32 | i; // modify the last
      }
      if (imax > gmax) {
         gmax = imax;
         te = i; // te is the end position on the target
         for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
            _mm256_store_si256(Hmax + j, _mm256_load_si256(H1 + j));
         if (gmax + q->shift >= 255 || gmax >= endsc)
            break;
      }
      S = H1;
      H1 = H0;
      H0 = S; // swap H0 and H1
   }
   r.score = gmax + q->shift < 255 ? gmax : 255;
   r.te = te;
   if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
      int max = -1, tmp, low, high, qlen = slen * 32;
      uint8_t *t = (uint8_t*) Hmax;
      for (i = 0; i < qlen; ++i, ++t)
         if ((int) *t > max)
            max = *t, r.qe = i / 32 + i % 32 * slen;
         else if ((int) *t == max && (tmp = i / 32 + i % 32 * slen) < r.qe)
            r.qe = tmp;
      //printf("%d,%d\n", max, gmax);
      if (b) {
         i = (r.score + q->max - 1) / q->max;
         low = te - i;
         high = te + i;
         for (i = 0; i < n_b; ++i) {
            int e = (int32_t) b[i];
            if ((e < low || e > high) && (int) (b[i] >> 32) > r.score2)
               r.score2 = b[i] >> 32, r.te2 = e;
         }
      }
   }
   free(b);

   return r;
}

kswr_t ksw_i16(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
   int slen, i,segNum;
  
   __m256i ninf, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E;
   kswr_t r;


   // initialization
   r = g_defr;
   ninf = _mm256_set1_epi16(NEG_INF);
   oe_del = _mm256_set1_epi16(-_o_del - _e_del);
   e_del = _mm256_set1_epi16(-_e_del);
   oe_ins = _mm256_set1_epi16(-_o_ins - _e_ins);
   e_ins = _mm256_set1_epi16(-_e_ins);
   H0 = q->H0;
   H1 = q->H1;
   E = q->E;
   slen = q->slen;
//-NW vars---------------
   int16_t* boundary =(int16_t*)malloc(sizeof(int16_t)*tlen+1);
/* initialize H and E */
   {
        for (i=0; LIKELY(i<slen); ++i) {
            __m256i_16_t h;
            __m256i_16_t e;
            for (segNum=0; LIKELY(segNum<16); ++segNum) {
                int64_t tmp = -_o_del -_e_del*(16*slen+i);
                h.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
                tmp = tmp - _o_del;
                e.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
            }
            _mm256_store_si256(H0+i, h.m);
            _mm256_store_si256(E+i, e.m);
            
        }
    }

    {
        boundary[0] = 0;
        for (i=1; LIKELY(i<=tlen); ++i) {
            int64_t tmp =-_o_del -_e_del*(i-1);
            boundary[i] = tmp < INT16_MIN ? INT16_MIN : tmp;
        }
    }
// the core loop
   for (i = 0; i < tlen; ++i) {
      int j, k;       
      __m256i e, t, h, f = ninf,*S = q->qp + target[i] * slen; // s is the 1st score vector
         
	/* load final segment of pvHStore and shift left by 2 bytes */     
	 h = _mm256_load_si256(H0 + slen -1); // h={2,5,8,11,14,17,-1,-1} in the above example
         h = _mm256_slli_si256_rpl(h, 2);

	/* insert upper boundary condition */
	 h = _mm256_insert_epi16_rpl(h, boundary[i], 0);

	/* inner loop to process the query sequence */
      for (j = 0; LIKELY(j < slen); ++j) {
         h = _mm256_adds_epi16(h, _mm256_load_si256(S + j));
         e = _mm256_load_si256(E + j);

	 /* Get max from vH, vE and vF. */
         h = _mm256_max_epi16(h, e);
         h = _mm256_max_epi16(h, f);

	 /* Save vH values. */
         _mm256_store_si256(H1 + j, h);

	 /* Update vE value. */
	 e = _mm256_add_epi16(e, e_del);
         t= _mm256_add_epi16(h, oe_del);
         e = _mm256_max_epi16(e, t);
         
        // e = _mm256_max_epi16(e, h);
         _mm256_store_si256(E + j, e);

	 /* Update vF value. */
         f = _mm256_add_epi16(f, e_ins);
	 t= _mm256_add_epi16(h, oe_ins);
         f = _mm256_max_epi16(f, t);
         

	   /* Load the next vH. */
         h = _mm256_load_si256(H0 + j);
      }
	
      for (k = 0; LIKELY(k < 16); ++k) {
        int64_t tmp = boundary[i+1]-_o_del;
        int16_t tmp2 = tmp < INT16_MIN ? INT16_MIN : tmp;
	 f = _mm256_slli_si256_rpl(f, 2);
         f = _mm256_insert_epi16_rpl(f, tmp2, 0);
         for (j = 0; LIKELY(j < slen); ++j) {
            h = _mm256_load_si256(H1 + j);
            h = _mm256_max_epi16(h, f);	     
            _mm256_store_si256(H1 + j, h);

            h = _mm256_add_epi16(h, oe_ins);	
            f = _mm256_add_epi16(f, e_ins);   
            if (UNLIKELY(!_mm256_movemask_epi8(_mm256_cmpgt_epi16(f, h))))
               goto end_loop8; 
         }
      }
end_loop8:

      S = H1;
      H1 = H0;
      H0 = S; // swap H0 and H1
}    
     

    /* extract last value from the last column */
       {
        __m256i vH = _mm256_load_si256(H0+slen-1);
            vH = _mm256_slli_si256_rpl (vH, 2);
            r.score=(int16_t) _mm256_extract_epi16_rpl (vH, 15)-1;
    }
   return r;
}


