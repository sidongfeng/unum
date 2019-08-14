/*
 * Copyright (c) 2016, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by
 * G. Scott Lloyd, lloyd23@llnl.gov. LLNL-CODE-704762. All rights reserved.
 * 
 * This file is part of Unum. For details, see
 * http://github.com/LLNL/unum
 * 
 * Please also read COPYING � Our Notice and GNU Lesser General Public 
 * License. 
 * 
 * This program is free software; you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License (as published by the 
 * Free Software Foundation) version 2.1 dated February 1999. 
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and 
 * conditions of the GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with this program; if not, write to the Free Software Foundation, 
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
 */

#ifndef GBND_H_
#define GBND_H_

#include "glayer.h" /* gbnd_s, end_t */

#if defined(__cplusplus)
extern "C" {
#endif

int ltgQ(const gbnd_s *g, const gbnd_s *h);
int gtgQ(const gbnd_s *g, const gbnd_s *h);
int neqgQ(const gbnd_s *g, const gbnd_s *h);
int samegQ(const gbnd_s *g, const gbnd_s *h);
int cmpgQ(const gbnd_s *g, end_t ge, const gbnd_s *h, end_t he);

int spanszerogQ(const gbnd_s *g);
void intersectgQ(gbnd_s *a, const gbnd_s *g, const gbnd_s *h);

void plusg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);
void minusg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);
void timesg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);
void divideg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);

void powg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);
void squareg(gbnd_s *a, const gbnd_s *g);
void sqrtg(gbnd_s *a, const gbnd_s *g);
void negateg(gbnd_s *a, const gbnd_s *g);
void absg(gbnd_s *a, const gbnd_s *g);
void expg(gbnd_s *a, const gbnd_s *g);
void logg(gbnd_s *a, const gbnd_s *g);
void cosg(gbnd_s *a, const gbnd_s *g);
void sing(gbnd_s *a, const gbnd_s *g);
void tang(gbnd_s *a, const gbnd_s *g);
void cotg(gbnd_s *a, const gbnd_s *g);

void ming(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);
void maxg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y);

int cliplg(gbnd_s *a, const gbnd_s *g, const gbnd_s *h);
int cliphg(gbnd_s *a, const gbnd_s *g, const gbnd_s *h);

// my own work
void mpf_log(mpf_t a,const mpf_t b, unsigned int prec);
void mpf_exp(mpf_t a,const mpf_t b, unsigned int prec);
void mpf_pow(mpf_t min, mpf_t max, const mpf_t x,const mpf_t y, unsigned int prec);
void mpf_pow_approx(mpf_t pow,const mpf_t x,const mpf_t y, unsigned int prec);
void mpf_root(mpf_t a,const mpf_t b, unsigned int root, unsigned int prec);
void mpf_mul_si(mpf_t result, const mpf_t op, signed int si);
int mpf_is_even(const mpf_t a);
void print_mpf(const mpf_t a);

void Taylor_Series_Log (mpf_t lg, mpf_t z, unsigned int prec);
void Taylor_Series_Exp(mpf_t ex, const mpf_t z, unsigned int prec);
void smallest_term (mpf_t min, int prec);
void factorial (mpf_t inv, unsigned int n);
unsigned int find_fact_digits(const int n);
unsigned int find_pow_digits(const mpf_t x, const mpf_t y,unsigned int prec);
void power(gnum_s *a, gnum_s *b, const gnum_s *x, const gnum_s *y, int *nan);
int exact(const gnum_s *x,const gnum_s *y);

void combinegb(gbnd_s *interval, const gbnd_s *x, const gbnd_s *y);
void ulp(mpf_t small, int prec);
void concat(char* result, const char *s1, const char *s2);
void substring(char *res, char *str, size_t begin, size_t end);
void dec2fra(fraction *a,const mpf_t b);

int nneqgQ(const gbnd_s *g, const gbnd_s *h);

#if defined (__cplusplus)
}
#endif

#endif /* GBND_H_ */
