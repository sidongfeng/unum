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

#include "gbnd.h"
#include "hlayer.h"

/* adding */
#include "stdio.h" 
#include "string.h" 
#include <stdlib.h>
#include "uenv.h"
#include <math.h>

/* scratchpad */

/* Test if interval g is strictly less than interval h. */

int ltgQ(const gbnd_s *g, const gbnd_s *h)
{
	return !(g->nan || h->nan) && cmp_gn(&g->r, RE, &h->l, LE) < 0;
}

/* Test if interval g is strictly greater than interval h. */

int gtgQ(const gbnd_s *g, const gbnd_s *h)
{
	return !(g->nan || h->nan) && cmp_gn(&g->l, LE, &h->r, RE) > 0;
}

/* Test if interval g is nowhere equal to interval h. */

int neqgQ(const gbnd_s *g, const gbnd_s *h)
{
	return !(g->nan || h->nan) && (ltgQ(g, h) || gtgQ(g, h));
}

/* Test if interval g is identical to interval h. */

int samegQ(const gbnd_s *g, const gbnd_s *h)
{
	return (g->nan && h->nan) ||
		(
			g->nan == h->nan &&
			mpf_cmp(g->l.f, h->l.f) == 0 &&
			g->l.inf == h->l.inf &&
			g->l.open == h->l.open &&
			mpf_cmp(g->r.f, h->r.f) == 0 &&
			g->r.inf == h->r.inf &&
			g->r.open == h->r.open
		);
}

/* Compare end points from intervals g and h. */

int cmpgQ(const gbnd_s *g, end_t ge, const gbnd_s *h, end_t he)
{
	if (g->nan || h->nan) return 0; /* FIXME: what about NaNs? */
	else return cmp_gn(ge==LE ? &g->l : &g->r, ge, he==LE ? &h->l : &h->r, he);
}

/* Test if interval g spans zero. */

int spanszerogQ(const gbnd_s *g)
{
	return ((!mpf_sgn(g->l.f) && !g->l.open) || (!mpf_sgn(g->r.f) && !g->r.open)) ||
	        (mpf_sgn(g->l.f) < 0 && mpf_sgn(g->r.f) > 0);
}

/* Add two general intervals. */

void plusg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y)
{
	/* If any value is NaN, the result is also NaN. */
	if (x->nan || y->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,plusg,y,gb);
		return;
	}
	a->nan = 0;

	/* Compute left endpoint: */
	/* Cases involving exact infinity or -infinity: */
	if (x->l.inf && mpf_sgn(x->l.f) < 0 && !x->l.open) {
		if (y->l.inf && mpf_sgn(y->l.f) > 0 && !y->l.open) {a->nan = 1;}
		else {mpf_set_si(a->l.f, -1); a->l.inf = 1; a->l.open = 0;}
	} else if (y->l.inf && mpf_sgn(y->l.f) < 0 && !y->l.open) {
		if (x->l.inf && mpf_sgn(x->l.f) > 0 && !x->l.open) {a->nan = 1;}
		else {mpf_set_si(a->l.f, -1); a->l.inf = 1; a->l.open = 0;}
	} else if (
		(x->l.inf && mpf_sgn(x->l.f) > 0 && !x->l.open) ||
		(y->l.inf && mpf_sgn(y->l.f) > 0 && !y->l.open)
		) {
		mpf_set_si(a->l.f, 1); a->l.inf = 1; a->l.open = 0;
	} else if (x->l.inf && mpf_sgn(x->l.f) < 0) {
		if (y->l.inf && mpf_sgn(y->l.f) > 0 && !y->l.open)
			{mpf_set_si(a->l.f, 1); a->l.inf = 1; a->l.open = 0;}
		else
			{mpf_set_si(a->l.f, -1); a->l.inf = 1; a->l.open = 1;}
	} else if (y->l.inf && mpf_sgn(y->l.f) < 0) {
		if (x->l.inf && mpf_sgn(x->l.f) > 0 && !x->l.open)
			{mpf_set_si(a->l.f, 1); a->l.inf = 1; a->l.open = 0;}
		else
			{mpf_set_si(a->l.f, -1); a->l.inf = 1; a->l.open = 1;}
	} else {
		/* What's left is the arithmetic case, done with extended precision. */
		mpf_add(a->l.f, x->l.f, y->l.f);
		a->l.inf = 0;
		a->l.open = x->l.open || y->l.open;
		/* NOTE: this range check is not in the Mathematica notebook 6-24-2015 */
#if 0
		{
			mpf_t tmpf;
			mpf_init2(tmpf, PBITS);
			mpf_neg(tmpf, maxreal);
			if (mpf_cmp(a->l.f, tmpf) < 0) {
				mpf_set_si(a->l.f, -1); a->l.inf = 1; a->l.open = 1;
			} else if (mpf_cmp(a->l.f, maxreal) > 0) {
				mpf_set(a->l.f, maxreal); a->l.inf = 0; a->l.open = 1;
			}
			mpf_clear(tmpf);
		}
#endif
	}

	/* Compute right endpoint, using similar logic: */
	/* Cases involving exact infinity or -infinity: */
	if (x->r.inf && mpf_sgn(x->r.f) < 0 && !x->r.open) {
		if (y->r.inf && mpf_sgn(y->r.f) > 0 && !y->r.open) {a->nan = 1;}
		else {mpf_set_si(a->r.f, -1); a->r.inf = 1; a->r.open = 0;}
	} else if (y->r.inf && mpf_sgn(y->r.f) < 0 && !y->r.open) {
		if (x->r.inf && mpf_sgn(x->r.f) > 0 && !x->r.open) {a->nan = 1;}
		else {mpf_set_si(a->r.f, -1); a->r.inf = 1; a->r.open = 0;}
	} else if (
		(x->r.inf && mpf_sgn(x->r.f) > 0 && !x->r.open) ||
		(y->r.inf && mpf_sgn(y->r.f) > 0 && !y->r.open)
		) {
		mpf_set_si(a->r.f, 1); a->r.inf = 1; a->r.open = 0;
	} else if (x->r.inf && mpf_sgn(x->r.f) > 0) {
		if (y->r.inf && mpf_sgn(y->r.f) < 0 && !y->r.open)
			{mpf_set_si(a->r.f, -1); a->r.inf = 1; a->r.open = 0;}
		else
			{mpf_set_si(a->r.f, 1); a->r.inf = 1; a->r.open = 1;}
	} else if (y->r.inf && mpf_sgn(y->r.f) > 0) {
		if (x->r.inf && mpf_sgn(x->r.f) < 0 && !x->r.open)
			{mpf_set_si(a->r.f, -1); a->r.inf = 1; a->r.open = 0;}
		else
			{mpf_set_si(a->r.f, 1); a->r.inf = 1; a->r.open = 1;}
	} else {
		/* What's left is the arithmetic case, done with extended precision. */
		mpf_add(a->r.f, x->r.f, y->r.f);
		a->r.inf = 0;
		a->r.open = x->r.open || y->r.open;
		/* TODO: this range check is not in the Mathematica notebook 6-24-2015 */
#if 0
		{
			mpf_t tmpf;
			mpf_init2(tmpf, PBITS);
			mpf_neg(tmpf, maxreal);
			if (mpf_cmp(a->r.f, tmpf) < 0) {
				mpf_set(a->r.f, tmpf); a->r.inf = 0; a->r.open = 1;
			} else if (mpf_cmp(a->r.f, maxreal) > 0) {
				mpf_set_si(a->r.f, 1); a->r.inf = 1; a->r.open = 1;
			}
			mpf_clear(tmpf);
		}
#endif
	}

	if (a->nan) {
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,plusg,y,gb);
	}
}

/* Subtraction in the g-layer. */

void minusg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y)
{
	gbnd_s gb;

	gbnd_init(&gb);
	negateg(&gb, y);
	plusg(a, x, &gb);
	gbnd_clear(&gb);
}

/* The "left" multiplication table for general intervals. */

static int timesposleft(gnum_s *a, const gnum_s *x, const gnum_s *y)
{
	int nan = 0;

	if (!mpf_sgn(x->f) && !x->open) {
		if (y->inf && !y->open) {nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
		else {mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;}
	} else if (!mpf_sgn(y->f) && !y->open) {
		if (x->inf && !x->open) {nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
		else {mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;}
	} else if (!mpf_sgn(x->f) && x->open) {
		if (y->inf && !y->open) {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;}
		else {mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
	} else if (!mpf_sgn(y->f) && y->open) {
		if (x->inf && !x->open) {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;}
		else {mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
	} else if ((x->inf && !x->open) || (y->inf && !y->open)) {
		mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;
	} else {
		mpf_mul(a->f, x->f, y->f);
		a->inf = 0;
		a->open = x->open || y->open;
	}
	return nan;
}

/* The "right" multiplication table for general intervals. */

static int timesposright(gnum_s *a, const gnum_s *x, const gnum_s *y)
{
	int nan = 0;

	if (x->inf && !x->open) {
		if (!mpf_sgn(y->f) && !y->open) {nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
		else {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;}
	} else if (y->inf && !y->open) {
		if (!mpf_sgn(x->f) && !x->open) {nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
		else {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;}
	} else if (x->inf && x->open) {
		if (!mpf_sgn(y->f) && !y->open) {mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;}
		else {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 1;}
	} else if (y->inf && y->open) {
		if (!mpf_sgn(x->f) && !x->open) {mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;}
		else {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 1;}
	} else if ((!mpf_sgn(x->f) && !x->open) || (!mpf_sgn(y->f) && !y->open)) {
		mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;
	} else {
		mpf_mul(a->f, x->f, y->f);
		a->inf = 0;
		a->open = x->open || y->open;
	}
	return nan;
}

/* Multiplication in the g-layer. */

void timesg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y)
{
	gnum_s lcan, rcan;
	gnum_s agn, xgn, ygn;

	/* If any value is NaN, the result is also NaN. */
	if (x->nan || y->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,timesg,y,gb);
		return;
	}
	a->nan = 0;

	gnum_init(&lcan);
	gnum_init(&rcan);
	gnum_init(&agn);
	gnum_init(&xgn);
	gnum_init(&ygn);

	/* Lower left corner is in upper right quadrant, facing uphill: */
	if (mpf_sgn(x->l.f) >= 0 && mpf_sgn(y->l.f) >= 0) {
		a->nan |= timesposleft(&lcan, &x->l, &y->l);
	} else {mpf_set_si(lcan.f, 1); lcan.inf = 1; lcan.open = 0;}
	/* Upper right corner is in lower left quadrant, facing uphill: */
	if ((mpf_sgn(x->r.f) < 0 || ((mpf_sgn(x->r.f) == 0) && x->r.open)) &&
	    (mpf_sgn(y->r.f) < 0 || ((mpf_sgn(y->r.f) == 0) && y->r.open))) {
		neg_gn(&xgn, &x->r);
		neg_gn(&ygn, &y->r);
		a->nan |= timesposleft(&agn, &xgn, &ygn);
		if (!a->nan && cmp_gn(&agn, LE, &lcan, LE) < 0) {
			mpf_set(lcan.f, agn.f);
			lcan.inf = agn.inf;
			lcan.open = agn.open;
		}
	}
	/* Upper left corner is in upper left quadrant, facing uphill: */
	if ((mpf_sgn(x->l.f) < 0 || ((mpf_sgn(x->l.f) == 0) && !x->l.open)) &&
	    (mpf_sgn(y->r.f) > 0 || ((mpf_sgn(y->r.f) == 0) && !y->r.open))) {
		neg_gn(&xgn, &x->l);
		a->nan |= timesposright(&agn, &xgn, &y->r);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, LE, &lcan, LE) < 0) {
			mpf_set(lcan.f, agn.f);
			lcan.inf = agn.inf;
			lcan.open = agn.open;
		}
	}
	/* Lower right corner is in lower right quadrant, facing uphill: */
	if ((mpf_sgn(x->r.f) > 0 || ((mpf_sgn(x->r.f) == 0) && !x->r.open)) &&
	    (mpf_sgn(y->l.f) < 0 || ((mpf_sgn(y->l.f) == 0) && !y->l.open))) {
		neg_gn(&ygn, &y->l);
		a->nan |= timesposright(&agn, &x->r, &ygn);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, LE, &lcan, LE) < 0) {
			mpf_set(lcan.f, agn.f);
			lcan.inf = agn.inf;
			lcan.open = agn.open;
		}
	}

	/* Upper right corner is in upper right quadrant, facing downhill: */
	if ((mpf_sgn(x->r.f) > 0 || ((mpf_sgn(x->r.f) == 0) && !x->r.open)) &&
	    (mpf_sgn(y->r.f) > 0 || ((mpf_sgn(y->r.f) == 0) && !y->r.open))) {
		a->nan |= timesposright(&rcan, &x->r, &y->r);
	} else {mpf_set_si(rcan.f, -1); rcan.inf = 1; rcan.open = 0;}
	/* Lower left corner is in lower left quadrant, facing downhill: */
	if ((mpf_sgn(x->l.f) < 0 || ((mpf_sgn(x->l.f) == 0) && !x->l.open)) &&
	    (mpf_sgn(y->l.f) < 0 || ((mpf_sgn(y->l.f) == 0) && !y->l.open))) {
		neg_gn(&xgn, &x->l);
		neg_gn(&ygn, &y->l);
		a->nan |= timesposright(&agn, &xgn, &ygn);
		if (!a->nan && cmp_gn(&agn, RE, &rcan, RE) > 0) {
			mpf_set(rcan.f, agn.f);
			rcan.inf = agn.inf;
			rcan.open = agn.open;
		}
	}
	/* Lower right corner is in upper left quadrant, facing downhill: */
	if ((mpf_sgn(x->r.f) < 0 || ((mpf_sgn(x->r.f) == 0) && x->r.open)) &&
	    (mpf_sgn(y->l.f) >= 0)) {
		neg_gn(&xgn, &x->r);
		a->nan |= timesposleft(&agn, &xgn, &y->l);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, RE, &rcan, RE) > 0) {
			mpf_set(rcan.f, agn.f);
			rcan.inf = agn.inf;
			rcan.open = agn.open;
		}
	}
	/* Upper left corner is in lower right quadrant, facing downhill: */
	if ((mpf_sgn(x->l.f) >= 0) &&
	    (mpf_sgn(y->r.f) < 0 || ((mpf_sgn(y->r.f) == 0) && y->r.open))) {
		neg_gn(&ygn, &y->r);
		a->nan |= timesposleft(&agn, &x->l, &ygn);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, RE, &rcan, RE) > 0) {
			mpf_set(rcan.f, agn.f);
			rcan.inf = agn.inf;
			rcan.open = agn.open;
		}
	}

	if (a->nan) {
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,timesg,y,gb);
	}	else {
		mpf_set(a->l.f, lcan.f);
		a->l.inf = lcan.inf;
		a->l.open = lcan.open;
		mpf_set(a->r.f, rcan.f);
		a->r.inf = rcan.inf;
		a->r.open = rcan.open;
	}
	gnum_clear(&ygn);
	gnum_clear(&xgn);
	gnum_clear(&agn);
	gnum_clear(&rcan);
	gnum_clear(&lcan);
}

/* The "left" division table for general intervals. */

static int divideposleft(gnum_s *a, const gnum_s *x, const gnum_s *y)
{
	int nan = 0;

	if (!mpf_sgn(y->f) && !y->open) {
		nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;
	} else if (x->inf && !x->open) {
		if (y->inf && !y->open) {nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
		else {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;}
	} else if ((!mpf_sgn(x->f) && !x->open) || (y->inf && !y->open)) {
		mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;
	} else if ((!mpf_sgn(x->f) && x->open) || (y->inf && y->open)) {
		mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;
	} else {
		mpf_div(a->f, x->f, y->f);
		a->inf = 0;
		a->open = x->open || y->open;
	}
	return nan;
}

/* The "right" division table for general intervals. */

static int divideposright(gnum_s *a, const gnum_s *x, const gnum_s *y)
{
	int nan = 0;

	if (!mpf_sgn(y->f) && !y->open) {
		nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;
	} else if (x->inf && !x->open) {
		if (y->inf && !y->open) {nan = 1; mpf_set_ui(a->f, 0); a->inf = 0; a->open = 1;}
		else {mpf_set_ui(a->f, 1); a->inf = 1; a->open = 0;}
	} else if ((!mpf_sgn(x->f) && !x->open) || (y->inf && !y->open)) {
		mpf_set_ui(a->f, 0); a->inf = 0; a->open = 0;
	} else if ((x->inf && x->open) || (!mpf_sgn(y->f) && y->open)) {
		mpf_set_ui(a->f, 1); a->inf = 1; a->open = 1;
	} else {
		mpf_div(a->f, x->f, y->f);
		a->inf = 0;
		a->open = x->open || y->open;
	}
	return nan;
}

/* Division in the g-layer. */

void divideg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y)
{
	gnum_s lcan, rcan;
	gnum_s agn, xgn, ygn;

	/* If any value is NaN, or denominator contains 0, the result is also NaN. */
	if (x->nan || y->nan ||
	    ((mpf_sgn(y->l.f) < 0 || ((mpf_sgn(y->l.f) == 0) && !y->l.open)) &&
	     (mpf_sgn(y->r.f) > 0 || ((mpf_sgn(y->r.f) == 0) && !y->r.open)))) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,divideg,y,gb);
		return;
	}
	a->nan = 0;

	gnum_init(&lcan);
	gnum_init(&rcan);
	gnum_init(&agn);
	gnum_init(&xgn);
	gnum_init(&ygn);

	/* FIXME: should not need a test for (y == 0(closed)), follows prototype */
	/* Upper left corner is in upper right quadrant, facing uphill: */
	if (mpf_sgn(x->l.f) >= 0 &&
	    (mpf_sgn(y->r.f) > 0 || ((mpf_sgn(y->r.f) == 0) && !y->r.open))) {
		a->nan |= divideposleft(&lcan, &x->l, &y->r);
	} else {mpf_set_si(lcan.f, 1); lcan.inf = 1; lcan.open = 0;}
	/* Lower right corner is in lower left quadrant, facing uphill: */
	if ((mpf_sgn(x->r.f) < 0 || ((mpf_sgn(x->r.f) == 0) &&  x->r.open)) &&
	    (mpf_sgn(y->l.f) < 0 || ((mpf_sgn(y->l.f) == 0) && !y->l.open))) {
		neg_gn(&xgn, &x->r);
		neg_gn(&ygn, &y->l);
		a->nan |= divideposleft(&agn, &xgn, &ygn);
		if (!a->nan && cmp_gn(&agn, LE, &lcan, LE) < 0) {
			mpf_set(lcan.f, agn.f);
			lcan.inf = agn.inf;
			lcan.open = agn.open;
		}
	}
	/* Lower left corner is in upper left quadrant, facing uphill: */
	if ((mpf_sgn(x->l.f) < 0 || ((mpf_sgn(x->l.f) == 0) && !x->l.open)) &&
	    mpf_sgn(y->l.f) >= 0) {
		neg_gn(&xgn, &x->l);
		a->nan |= divideposright(&agn, &xgn, &y->l);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, LE, &lcan, LE) < 0) {
			mpf_set(lcan.f, agn.f);
			lcan.inf = agn.inf;
			lcan.open = agn.open;
		}
	}
	/* Upper right corner is in lower right quadrant, facing uphill: */
	if ((mpf_sgn(x->r.f) > 0 || ((mpf_sgn(x->r.f) == 0) && !x->r.open)) &&
	    (mpf_sgn(y->r.f) < 0 || ((mpf_sgn(y->r.f) == 0) &&  y->r.open))) {
		neg_gn(&ygn, &y->r);
		a->nan |= divideposright(&agn, &x->r, &ygn);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, LE, &lcan, LE) < 0) {
			mpf_set(lcan.f, agn.f);
			lcan.inf = agn.inf;
			lcan.open = agn.open;
		}
	}

	/* Lower right corner is in upper right quadrant, facing downhill: */
	if ((mpf_sgn(x->r.f) > 0 || ((mpf_sgn(x->r.f) == 0) && !x->r.open)) &&
	    mpf_sgn(y->l.f) >= 0) {
		a->nan |= divideposright(&rcan, &x->r, &y->l);
	} else {mpf_set_si(rcan.f, -1); rcan.inf = 1; rcan.open = 0;}
	/* Upper left corner is in lower left quadrant, facing downhill: */
	if ((mpf_sgn(x->l.f) < 0 || ((mpf_sgn(x->l.f) == 0) && !x->l.open)) &&
	    (mpf_sgn(y->r.f) < 0 || ((mpf_sgn(y->r.f) == 0) &&  y->r.open))) {
		neg_gn(&xgn, &x->l);
		neg_gn(&ygn, &y->r);
		a->nan |= divideposright(&agn, &xgn, &ygn);
		if (!a->nan && cmp_gn(&agn, RE, &rcan, RE) > 0) {
			mpf_set(rcan.f, agn.f);
			rcan.inf = agn.inf;
			rcan.open = agn.open;
		}
	}
	/* Upper right corner is in upper left quadrant, facing downhill: */
	if ((mpf_sgn(x->r.f) < 0 || ((mpf_sgn(x->r.f) == 0) &&  x->r.open)) &&
	    (mpf_sgn(y->r.f) > 0 || ((mpf_sgn(y->r.f) == 0) && !y->r.open))) {
		neg_gn(&xgn, &x->r);
		a->nan |= divideposleft(&agn, &xgn, &y->r);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, RE, &rcan, RE) > 0) {
			mpf_set(rcan.f, agn.f);
			rcan.inf = agn.inf;
			rcan.open = agn.open;
		}
	}
	/* Lower left corner is in lower right quadrant, facing downhill: */
	if (mpf_sgn(x->l.f) >= 0 &&
	    (mpf_sgn(y->l.f) < 0 || ((mpf_sgn(y->l.f) == 0) && !y->l.open))) {
		neg_gn(&ygn, &y->l);
		a->nan |= divideposleft(&agn, &x->l, &ygn);
		neg_gn(&agn, &agn);
		if (!a->nan && cmp_gn(&agn, RE, &rcan, RE) > 0) {
			mpf_set(rcan.f, agn.f);
			rcan.inf = agn.inf;
			rcan.open = agn.open;
		}
	}

	if (a->nan) {
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,divideg,y,gb);
	}	else {
		mpf_set(a->l.f, lcan.f);
		a->l.inf = lcan.inf;
		a->l.open = lcan.open;
		mpf_set(a->r.f, rcan.f);
		a->r.inf = rcan.inf;
		a->r.open = rcan.open;
	}
	gnum_clear(&ygn);
	gnum_clear(&xgn);
	gnum_clear(&agn);
	gnum_clear(&rcan);
	gnum_clear(&lcan);
}

/* Square in the g-layer. */

void squareg(gbnd_s *a, const gbnd_s *g)
{
	gnum_s t1, t2;
	gnum_s *aL, *aR;

	if (g->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,squareg,g,gb);
		return;
	}
	a->nan = 0;

	gnum_init(&t1);
	gnum_init(&t2);

	mpf_pow_ui(t1.f, g->l.f, 2);
	t1.inf = g->l.inf;
	t1.open = g->l.open;

	mpf_pow_ui(t2.f, g->r.f, 2);
	t2.inf = g->r.inf;
	t2.open = g->r.open;

	if (cmp_gn(&t1, RE, &t2, RE) > 0) {
		aL = &t2; aR = &t1;
	} else {
		aL = &t1; aR = &t2;
	}

	/* See if 0 is in the range */
	if (spanszerogQ(g)) {
		mpf_set_si(a->l.f, 0);
		a->l.inf = 0;
		a->l.open = 0;
	} else {
		mpf_set(a->l.f, aL->f);
		a->l.inf = aL->inf;
		a->l.open = aL->open;
	}

	mpf_set(a->r.f, aR->f);
	a->r.inf = aR->inf;
	a->r.open = aR->open;

	gnum_clear(&t2);
	gnum_clear(&t1);
}

/* Square root in the g-layer. */

void sqrtg(gbnd_s *a, const gbnd_s *g)
{
	if (g->nan || mpf_sgn(g->l.f) < 0) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,sqrtg,g,gb);
		return;
	}
	a->nan = 0;

	mpf_sqrt(a->l.f, g->l.f);
	mpf_sqrt(a->r.f, g->r.f);
	a->l.inf = g->l.inf;
	a->r.inf = g->r.inf;
	a->l.open = g->l.open;
	a->r.open = g->r.open;
}

/* Negate a general interval. */

void negateg(gbnd_s *a, const gbnd_s *g)
{
	gnum_s tmp;
	const gnum_s *gl, *gr;

	/* Handle NaN input */
	if (g->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,negateg,g,gb);
		return;
	}
	a->nan = 0;

	/* Use temporary if input and output are the same */
	if (a == g) {
		gnum_init(&tmp);
		mpf_set(tmp.f, g->l.f);
		tmp.inf  = g->l.inf;
		tmp.open = g->l.open;
		gl = &tmp;
	} else {
		gl = &g->l;
	}
	gr = &g->r;

	/* Negate and reverse */
	mpf_neg(a->l.f, gr->f);
	a->l.inf  = gr->inf;
	a->l.open = gr->open;

	mpf_neg(a->r.f, gl->f);
	a->r.inf  = gl->inf;
	a->r.open = gl->open;

	if (a == g) gnum_clear(&tmp);
}

/* Absolute value in the g-layer. */

void absg(gbnd_s *a, const gbnd_s *g)
{
	gbnd_s tmp;

	if (g->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,absg,g,gb);
		return;
	}
	a->nan = 0;

	gbnd_init(&tmp);
	mpf_abs(tmp.l.f, g->l.f);
	tmp.l.inf = g->l.inf;
	tmp.l.open = g->l.open;
	mpf_abs(tmp.r.f, g->r.f);
	tmp.r.inf = g->r.inf;
	tmp.r.open = g->r.open;
	if (mpf_sgn(g->r.f) <= 0) {
		mpf_set(a->l.f, tmp.r.f);
		a->l.inf = tmp.r.inf;
		a->l.open = tmp.r.open;
		mpf_set(a->r.f, tmp.l.f);
		a->r.inf = tmp.l.inf;
		a->r.open = tmp.l.open;
	} else if (mpf_sgn(g->l.f) <= 0) {
		mpf_set_si(a->l.f, 0);
		a->l.inf = 0;
		a->l.open = 0;
		if (mpf_cmp(tmp.l.f, tmp.r.f) < 0) {
			mpf_set(a->r.f, tmp.r.f);
			a->r.inf = tmp.r.inf;
			a->r.open = tmp.r.open;
		} else if (mpf_cmp(tmp.l.f, tmp.r.f) > 0) {
			mpf_set(a->r.f, tmp.l.f);
			a->r.inf = tmp.l.inf;
			a->r.open = tmp.l.open;
		} else {
			mpf_set(a->r.f, tmp.r.f);
			a->r.inf = tmp.r.inf;
			a->r.open = tmp.l.open & tmp.r.open;
		}
	} else {
		mpf_set(a->l.f, tmp.l.f);
		a->l.inf = tmp.l.inf;
		a->l.open = tmp.l.open;
		mpf_set(a->r.f, tmp.r.f);
		a->r.inf = tmp.r.inf;
		a->r.open = tmp.r.open;
	}

	gbnd_clear(&tmp);
}

void ming(gbnd_s *a, const gbnd_s *x, const gbnd_s *y)
{
	/* If any value is NaN, the result is also NaN. */
	if (x->nan || y->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,ming,y,gb);
		return;
	}
	a->nan = 0;

	if (cmp_gn(&x->l, LE, &y->l, LE) < 0) {
		mpf_set(a->l.f, x->l.f);
		a->l.inf = x->l.inf;
		a->l.open = x->l.open;
	} else {
		mpf_set(a->l.f, y->l.f);
		a->l.inf = y->l.inf;
		a->l.open = y->l.open;
	}
	if (cmp_gn(&x->r, RE, &y->r, RE) < 0) {
		mpf_set(a->r.f, x->r.f);
		a->r.inf = x->r.inf;
		a->r.open = x->r.open;
	} else {
		mpf_set(a->r.f, y->r.f);
		a->r.inf = y->r.inf;
		a->r.open = y->r.open;
	}
}

void maxg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y)
{
	/* If any value is NaN, the result is also NaN. */
	if (x->nan || y->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,maxg,y,gb);
		return;
	}
	a->nan = 0;

	if (cmp_gn(&x->l, LE, &y->l, LE) > 0) {
		mpf_set(a->l.f, x->l.f);
		a->l.inf = x->l.inf;
		a->l.open = x->l.open;
	} else {
		mpf_set(a->l.f, y->l.f);
		a->l.inf = y->l.inf;
		a->l.open = y->l.open;
	}
	if (cmp_gn(&x->r, RE, &y->r, RE) > 0) {
		mpf_set(a->r.f, x->r.f);
		a->r.inf = x->r.inf;
		a->r.open = x->r.open;
	} else {
		mpf_set(a->r.f, y->r.f);
		a->r.inf = y->r.inf;
		a->r.open = y->r.open;
	}
}

int cliplg(gbnd_s *a, const gbnd_s *g, const gbnd_s *h)
{
	int res = 0;

	/* If any value is NaN, the result is also NaN. */
	if (g->nan || h->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,g,gb,cliplg,h,gb);
		return 0;
	}
	a->nan = 0;

	if (cmp_gn(&g->l, LE, &h->l, LE) < 0) {
		mpf_set(a->l.f, h->l.f);
		a->l.inf = h->l.inf;
		a->l.open = h->l.open;
		res = 1;
	} else {
		mpf_set(a->l.f, g->l.f);
		a->l.inf = g->l.inf;
		a->l.open = g->l.open;
	}
	if (cmp_gn(&g->r, RE, &h->r, RE) < 0) {
		mpf_set(a->r.f, h->r.f);
		a->r.inf = h->r.inf;
		a->r.open = h->r.open;
		res = 1;
	} else {
		mpf_set(a->r.f, g->r.f);
		a->r.inf = g->r.inf;
		a->r.open = g->r.open;
	}
	return res;
}

int cliphg(gbnd_s *a, const gbnd_s *g, const gbnd_s *h)
{
	int res = 0;

	/* If any value is NaN, the result is also NaN. */
	if (g->nan || h->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,g,gb,cliphg,h,gb);
		return 0;
	}
	a->nan = 0;

	if (cmp_gn(&g->l, LE, &h->l, LE) > 0) {
		mpf_set(a->l.f, h->l.f);
		a->l.inf = h->l.inf;
		a->l.open = h->l.open;
		res = 1;
	} else {
		mpf_set(a->l.f, g->l.f);
		a->l.inf = g->l.inf;
		a->l.open = g->l.open;
	}
	if (cmp_gn(&g->r, RE, &h->r, RE) > 0) {
		mpf_set(a->r.f, h->r.f);
		a->r.inf = h->r.inf;
		a->r.open = h->r.open;
		res = 1;
	} else {
		mpf_set(a->r.f, g->r.f);
		a->r.inf = g->r.inf;
		a->r.open = g->r.open;
	}
	return res;
}

// -----------------------------------------------
//        below is my own work
// -----------------------------------------------


/*  @ parameter     output(mpf_t), input(mpf_t), precision(unsigned int)
	Algorithm:	Find a equation log(x) = n*log(2)+log(b) where 0< b <= 2
				log(4.8) = 2*log(2) + log(1.2)
				calculate left part by brute force implementation,
						  right part by Taylor Expansion
*/
void mpf_log(mpf_t log, const mpf_t x, unsigned int prec){
	// error when x<0
	if (mpf_cmp_ui(x,0)<=0){
		printf("Error!\n");
		return;
	}
	
	// ln(1) = 0
	if (mpf_cmp_ui(x,1)==0){
		mpf_set_ui(log,0);
		return;
	}
	mp_bitcnt_t bits = prec * 3.322 + 4;
	mpf_t tmp,ln2;
	mpf_init2(tmp,bits);
	mpf_init2(ln2,bits);
	mpf_set(tmp,x);
	mpf_set_str(ln2,"0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481020570685733685520235758130557032",10);
	int count_n = 0;
	/* count how many ln(2)  */
	while (mpf_cmp_ui(tmp, 2) > 0){
		mpf_div_ui (tmp, tmp, 2);
		count_n += 1;
	}
	while (mpf_cmp_ui(tmp, 1) < 0){
		mpf_mul_ui (tmp, tmp, 2);
		count_n -= 1;
	}
	
	/* calculate other part by using Taylor Expansion */
	if (!mpf_cmp_ui(tmp,2)){
		mpf_set(log,ln2);
	}else if (mpf_cmp_ui(tmp,1)==0){
		mpf_set_ui(log,0);
	}else{
		mpf_ui_sub(tmp, 1, tmp);
		Taylor_Series_Log(log, tmp, prec+1);
		mpf_neg(log, log);
	}
	if (count_n){
		/* implement mpf_mul_si
		   add the n*log(2)+log(b) as the result
		*/
		mpf_mul_si(tmp,ln2,count_n);
		mpf_add (log,log,tmp);
	}

	mpf_clear (ln2);
	mpf_clear (tmp);
}

/* exp(real) * exp(decimal) */
void mpf_exp(mpf_t exp, const mpf_t x, unsigned int prec){
	/* exp(0) = 1 */
	if (!mpf_cmp_ui(x,0)){
		mpf_set_ui(exp,1);
		return;
	}
	mp_bitcnt_t bits = prec * 3.322+4;
	mpf_t tmp, e, u;
	mpf_init2 (tmp, bits);
	mpf_init2 (e, bits);
	mpf_init2 (u, bits);
	ulp(u,prec+1);
	mpf_set_str (e, "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290",10);
	
	/* find real and decimal */
	signed int real = mpf_get_si (x);
	if (real){
		if (real>0){
			mpf_sub_ui(tmp,x,real);
		}else{
			mpf_add_ui(tmp,x,-real);
		}
	}else{
		mpf_set(tmp,x);
	}

	Taylor_Series_Exp(exp, tmp, prec+1);
	
	if (real){
		if (real>0){
			mpf_pow_ui (tmp, e, real);
		}else{
			mpf_add(e,e,u);
			mpf_pow_ui (tmp, e, -real);
			mpf_ui_div(tmp,1,tmp);
		}
		mpf_mul (exp, exp, tmp);
	}

	mpf_clear (tmp);
	mpf_clear (e);
	mpf_clear(u);
}

/* calculate power into an interval without any loss in accuracy
	the result needs to be correct in any circumstance
*/
void mpf_pow(mpf_t min, mpf_t max, const mpf_t x,const mpf_t y, unsigned int prec){
	mp_bitcnt_t bits = prec * 3.322+4;
	mpf_t tmp,u;
	mpf_init2(tmp,bits);
	mpf_init2(u,bits);

	/* (anything except 0) ^ 0  is 1 */
	if (!mpf_cmp_ui(y,0) && mpf_cmp_ui(x,0)){mpf_set_ui(max,1);mpf_set_ui(min,1);return;}
	/* (anything) ^ 1  is itself */
	if (!mpf_cmp_ui(y,1)){mpf_set(max,x);mpf_set(min,x);return;}
	/* 1 ^ anything is 1 */
	if (!mpf_cmp_ui(x,1)){mpf_set_ui(max,1);mpf_set_ui(min,1);return;}
	/* 0 ^ +anything is 0
	       -anything and 0 is undefined
	*/
	if (!mpf_cmp_ui(x,0)){
		if (mpf_cmp_ui(y,0)>0){
			mpf_set_ui(max,0);
			mpf_set_ui(min,0);
			return;
		}else{
			//printf("Undefined");
			return;
		}
	}
	/* (-anything) ^ (non integer) is undefined */
	if (mpf_cmp_ui(x,0)<0 && !mpf_integer_p(y)){
		//printf("Undefined");
		return;
	}

	/* if the exponenet is integer, calculate by mpf_pow_ui */
	if (mpf_integer_p(y)){
		unsigned int i = mpf_get_ui(y);
		if (mpf_cmp_ui(y,0)>0){
			mpf_pow_ui(tmp,x,i);
		}else{
			mpf_pow_ui(tmp,x,i);
			mpf_ui_div(tmp,1,tmp);
		}
		mpf_set(max,tmp);
		mpf_set(min,tmp);
		return;
	}

	/* if the exponenet can be square rooted, calculate by mpf_sqrt */
	mpf_mul_ui(tmp,y,2);
	if (mpf_integer_p(tmp)){
		mpf_sqrt(min,x);
		mpf_floor(tmp,y);
		unsigned int i = mpf_get_ui(tmp);
		if (mpf_cmp_ui(y,0)>0){
			mpf_pow_ui(tmp,x,i);
		}else{
			mpf_pow_ui(tmp,x,i);
			mpf_ui_div(tmp,1,tmp);
		}
		mpf_mul(tmp,tmp,min);
		mpf_set(max,tmp);
		mpf_set(min,tmp);
		return;
	}

	ulp(u,prec);
	/* find the upper bound and lower bound */

	/* x^y = exp(y*log(x)) */
	mpf_set(tmp,x);
	mpf_abs(tmp,tmp);
	mpf_log(tmp,tmp,prec+1);
	mpf_mul(tmp,y,tmp);
	mpf_exp(tmp,tmp,prec+1);
	/* (-anything) ^ (integer) is (-1)^(integer)*(+anything^integer) */
	if (!mpf_is_even(y) && mpf_cmp_ui(x,0)<0) mpf_neg(tmp,tmp);
	mpf_set(min,tmp);
	mpf_set(max,tmp);
	
	/* x^y = exp(y*log(x+ulp)) */
	mpf_set(tmp,x);
	mpf_abs(tmp,tmp);
	mpf_add(tmp,tmp,u);
	mpf_log(tmp,tmp,prec+1);
	mpf_mul(tmp,y,tmp);
	mpf_exp(tmp,tmp,prec+1);
	if (!mpf_is_even(y) && mpf_cmp_ui(x,0)<0) mpf_neg(tmp,tmp);
	if (mpf_cmp(min,tmp)>0) mpf_set(min,tmp);
	if (mpf_cmp(max,tmp)<0) mpf_set(max,tmp);

	/* x^y = exp(y*log(x)+ulp) */
	mpf_set(tmp,x);
	mpf_abs(tmp,tmp);
	mpf_log(tmp,tmp,prec+1);
	mpf_mul(tmp,y,tmp);
	mpf_add(tmp,tmp,u);
	mpf_exp(tmp,tmp,prec+1);
	if (!mpf_is_even(y) && mpf_cmp_ui(x,0)<0) mpf_neg(tmp,tmp);
	if (mpf_cmp(min,tmp)>0) mpf_set(min,tmp);
	if (mpf_cmp(max,tmp)<0) mpf_set(max,tmp);

	/* x^y = exp(y*log(x+ulp)+ulp) */
	mpf_set(tmp,x);
	mpf_abs(tmp,tmp);
	mpf_add(tmp,tmp,u);
	mpf_log(tmp,tmp,prec+1);
	mpf_mul(tmp,y,tmp);
	mpf_add(tmp,tmp,u);
	mpf_exp(tmp,tmp,prec+1);
	if (!mpf_is_even(y) && mpf_cmp_ui(x,0)<0) mpf_neg(tmp,tmp);
	if (mpf_cmp(min,tmp)>0) mpf_set(min,tmp);
	if (mpf_cmp(max,tmp)<0) mpf_set(max,tmp);

	mpf_clear(tmp);
	mpf_clear(u);
}

/* calculate x^y by using the formulae
	x^y = exp(y*log(x))
	but it loss accuracy due to exponential function suffers from error magnification
*/
void mpf_pow_approx(mpf_t pow, const mpf_t x,const mpf_t y, unsigned int prec){
	mp_bitcnt_t bits = prec * 3.322;
	mpf_t tmp;
	mpf_init2(tmp,bits);

	// (anything except 0) ^ 0  is 1
	if (!mpf_cmp_ui(y,0) && mpf_cmp_ui(x,0)){mpf_set_ui(pow,1);return;}
	// (anything) ^ 1  is itself
	if (!mpf_cmp_ui(y,1)){mpf_set(pow,x);return;}
	//  1 ^ (anything) is 1
	if (!mpf_cmp_ui(x,1)){mpf_set_ui(pow,1);return;}
	// 0 ^  +anything is 0
	//      -anything and 0 is undefined
	if (!mpf_cmp_ui(x,0)){
		if (mpf_cmp_ui(y,0)>0){
			mpf_set_ui(pow,0);
			return;
		}else{
			printf("Undefined");
			return;
		}
	}
	// (-anything) ^ (non integer) is undefined
	if (mpf_cmp_ui(x,0)<0 && !mpf_integer_p(y)){
		printf("Undefined");
		return;
	}

	/* x^y = exp(y*log(x)) */
	mpf_set(tmp,x);
	mpf_abs(tmp,tmp);
	mpf_log(tmp,tmp,prec);
	mpf_mul(tmp,y,tmp);
	mpf_exp(tmp,tmp,prec);
	if (!mpf_is_even(y) && mpf_cmp_ui(x,0)<0) mpf_neg(tmp,tmp);
	mpf_set(pow,tmp);

	mpf_clear(tmp);
}

/* nth root of a number (with bugs)*/
void mpf_root(mpf_t a,const mpf_t b, unsigned int root, unsigned int prec){
	/* 1 root of any number is itself */
	if (root==1){
		mpf_set(a,b);
		return;
	}
	mp_bitcnt_t bits = prec * 3.322+100;
	mpf_t dec,tmp,lhs,rhs,b_root;
	mpz_t real;
	mpz_init2(real,bits);
	mpf_init2(dec,bits);
	mpf_init2(tmp,bits);

	/* convert mpf_t to a char list in decimal */
	/* 123.456 is converted to real = 123, dec = 456 */
	mpz_set_f(real,b);
	mpf_set_z(tmp,real);
	mpf_sub(dec,b,tmp);
	char real_s[256],dec_s[256];
	mp_exp_t dec_e;
	mpz_get_str (real_s,10,real);
	char dectmp[256];
	mpf_get_str (dectmp,&dec_e,10,0,dec);
	if (dec_e == 0){
		mpf_get_str (dec_s,&dec_e,10,0,dec);
	}else{
		char zero[-dec_e];
		for (int i=0; i<-dec_e;i++){ 
			zero[i] = '0'; 
		}
		concat(dec_s,zero,dectmp);
	}
	//printf("real: %s\n",real_s);
	//printf("dec: %s\n",dec_s);
	
	mpf_init2(lhs,bits);
	mpf_init2(rhs,bits);
	mpf_init2(b_root,bits);
	mpf_set_str(b_root,"10",10);
	mpf_pow_ui(b_root,b_root,root);

	int rptr = 0;
	int dptr = 0;
	mpf_t r,y,y_l,alp;
	mpf_init2(r,bits);
	mpf_init2(y,bits);
	mpf_init2(y_l,bits);
	mpf_init2(alp,bits);
	int beta = 0;
	mpf_set_ui(r,0);
	mpf_set_ui(y,0);
	mpf_set_ui(y_l,0);

	int c = 0;
	int end_exact = 0;

	char result_r[(int)ceil(strlen(real_s)/root)];
	char result_d[prec-(int)ceil(strlen(real_s)/root)];
	int result_r_ptr = 0;
	int result_d_ptr = 0;

	/* terminate until find exact result, or in the last expected decimal */
	while (!(!(c < prec) || end_exact)){
		/* check r, alpha, beta, y and y_l */
		/*
		printf("check: ");
		printf("r: ");print_mpf(r);
		printf("alpha: ");print_mpf(alp);
		printf("beta: %i\n",beta);
		printf("y: ");print_mpf(y);
		printf("y_l: ");print_mpf(y_l);
		printf("\n");
		*/
		if (rptr < strlen(real_s)){
			// find alpha
			if (strlen(real_s) % root != 0 && rptr == 0){
				char alpha[strlen(real_s) % root];
				substring(alpha,real_s,rptr,strlen(real_s) % root);
				mpf_set_str(alp,alpha,10);
				rptr = strlen(real_s) % root;
			}else{
				char alpha[root];
				substring(alpha,real_s,rptr,rptr+root);
				mpf_set_str(alp,alpha,10);
				rptr += root;
			}

			// find beta
			beta = 0;
			while (1){
				mpf_mul_ui(lhs,y,10);
				mpf_add_ui(lhs,lhs,beta);
				mpf_pow_ui(lhs,lhs,root);
				mpf_pow_ui(rhs,y,root);
				mpf_add(rhs,rhs,r);
				mpf_mul(rhs,rhs,b_root);
				mpf_add(rhs,rhs,alp);
				if (mpf_cmp(lhs,rhs) > 0) break;
				beta ++;
			}
			beta--;
			
			// find y_l
			mpf_set(y_l,y);

			// find y
			mpf_mul_ui(tmp,y,10);
			mpf_add_ui(y,tmp,beta);

			// find r
			mpf_mul(lhs,b_root,r);
			mpf_add(lhs,lhs,alp);
			mpf_pow_ui(tmp,y,root);
			mpf_sub(lhs,lhs,tmp);
			mpf_pow_ui(tmp,y_l,root);
			mpf_mul(tmp,tmp,b_root);
			mpf_add(r,lhs,tmp);

			result_r[result_r_ptr] = beta +'0';
			result_r_ptr++;
		}else{
			// find alpha
			if (dptr < strlen(dec_s)){
				if (dptr + root <= strlen(dec_s)){
					char alpha[root];
					substring(alpha,dec_s,dptr,dptr+root);
					mpf_set_str(alp,alpha,10);
				}else{
					char alpha[root];
					substring(alpha,dec_s,dptr,strlen(dec_s));
					mpf_set_str(alp,alpha,10);
					mpf_pow_ui(alp,alp,root-(strlen(dec_s)-dptr));
				}
				dptr += root;
			}else{
				mpf_set_ui(alp,0);
			}

			beta = 0;
			while (1){
				mpf_mul_ui(lhs,y,10);
				mpf_add_ui(lhs,lhs,beta);
				mpf_pow_ui(lhs,lhs,root);
				mpf_pow_ui(rhs,y,root);
				mpf_add(rhs,rhs,r);
				mpf_mul(rhs,rhs,b_root);
				mpf_add(rhs,rhs,alp);
				if (mpf_cmp(lhs,rhs) > 0) break;
				beta ++;
			}
			beta--;
			
			// find y_l
			mpf_set(y_l,y);

			// find y
			mpf_mul_ui(tmp,y,10);
			mpf_add_ui(y,tmp,beta);

			// find r
			mpf_mul(lhs,b_root,r);
			mpf_add(lhs,lhs,alp);
			mpf_pow_ui(tmp,y,root);
			mpf_sub(lhs,lhs,tmp);
			mpf_pow_ui(tmp,y_l,root);
			mpf_mul(tmp,tmp,b_root);
			mpf_add(r,lhs,tmp);

			result_d[result_d_ptr] = beta +'0';
			result_d_ptr++;
		}
		end_exact = (dptr >= strlen(dec_s)) && !mpf_cmp_ui(r,0);
		c += 1;	
	}
	/* return result in decimal format */
	if (strlen(result_d)==0){
		char res[strlen(result_r)];
		strcpy(res,result_r);
		mpf_set_str(a,res,10);
	}else{
		char res[strlen(result_r)+strlen(result_d)+1];
		strcpy(res,result_r);
		strcat(res,".");
		strcat(res,result_d);
		mpf_set_str(a,res,10);
	}

	mpz_clear(real);
	mpf_clear(dec);
	mpf_clear(tmp);
	mpf_clear(lhs);
	mpf_clear(rhs);
	mpf_clear(b_root);
	mpf_clear(r);
	mpf_clear(y);
	mpf_clear(y_l);
	mpf_clear(alp);
}

/* multiple mpf_t with a signed integer */
void mpf_mul_si(mpf_t result, const mpf_t op, signed int si){
	if (!si || !mpf_cmp_ui(op,0)){
		mpf_set_ui(result,0);
		return;
	} 
	mpf_mul_ui (result, op, abs(si));
	if (si<0) mpf_neg (result, result);
}

/* check the number is even, 1 indicates even */
int mpf_is_even(const mpf_t a){
	signed int i= mpf_get_si(a);
	return !(i % 2);
}

/* print out mpf value */
void print_mpf(const mpf_t a){
	char tmp[256];
	mp_exp_t tmp_exp;
	mpf_get_str(tmp,&tmp_exp,10,0,a);
	printf("mpf: %s with %li digits in real\n",tmp,tmp_exp);
}

/* log(a,b) = (log(a),log(b)) */
void logg(gbnd_s *a, const gbnd_s *g){

	gnum_s t1, t2;
	
	// if it is NaN, it returns NaN
	if (g->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,absg,g,gb);
		return;
	}
	// if any numbers is less than 0, it returns NaN
	if (mpf_cmp_ui(g->l.f,0)<0 || mpf_cmp_ui(g->r.f,0)<0){
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,absg,g,gb);
		return;
	}

	a->nan = 0;
	gnum_init(&t1);
	gnum_init(&t2);

	/* Evaluate left handside
		l = 0   ->  -inf
		  = inf ->  inf
		  = log(l)
	*/
	if (!mpf_cmp_ui(g->l.f,0)){
		mpf_set_si(t1.f,-1);
		t1.inf = 1;
		t1.open = g->l.open;
	}else if (g->l.inf){
		mpf_set_ui(t1.f,1);
		t1.inf = 1;
		t1.open = g->l.open;
	}else{
		mpf_log(t1.f, g->l.f,(PBITS+4)/3.322);
		t1.inf = g->l.inf;
		t1.open = g->l.open;
	}
	/* Evaluate right handside
		r = 0   ->  -inf
		  = inf ->  inf
		  = log(l)
	*/
	if (!mpf_cmp_ui(g->r.f,0)){
		mpf_set_si(t2.f,-1);
		t2.inf = 1;
		t2.open = g->r.open;
	}else if (g->r.inf){
		mpf_set_ui(t2.f,1);
		t2.inf = 1;
		t2.open = g->r.open;
	}else{
		mpf_log(t2.f, g->r.f,(PBITS+4)/3.322);
		t2.inf = g->r.inf;
		t2.open = g->r.open;
	}
	mpf_set(a->l.f, t1.f);
	a->l.inf = t1.inf;
	a->l.open = t1.open;
	mpf_set(a->r.f, t2.f);
	a->r.inf = t2.inf;
	a->r.open = t2.open;

	gnum_clear(&t2);
	gnum_clear(&t1);
	/* check gb */
	#if 0
	{
		printf("check gb: ");print_gb(a);printf("\n");
	}
	#endif
}

/* exp(a,b) = (exp(a),exp(b)) */
void expg(gbnd_s *a, const gbnd_s *g){
	gnum_s t1, t2;
	// NAN
	if (g->nan) {
		a->nan = 1;
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP1("NaN",a,gb,absg,g,gb);
		return;
	}

	a->nan = 0;
	gnum_init(&t1);
	gnum_init(&t2);
	/* calculate left endpoint */
	/* exp(x) = inf if x = inf
			  = 0   if x = -inf
			  = mpf_exp
	*/
	if (g->l.inf && mpf_sgn(g->l.f)>0){
		mpf_set_ui(a->l.f,1);
		a->l.inf = 1;
		a->l.open = g->l.open;
	}else if (g->l.inf && mpf_sgn(g->l.f)<0){
		mpf_set_ui(a->l.f,0);
		a->l.inf = 0;
		a->l.open = g->l.open;
	}else{
		mpf_exp(a->l.f, g->l.f,(PBITS+4)/3.322);
		a->l.inf = 0;
		a->l.open = g->l.open;
	}
	/* calculate right endpoint */
	if (g->r.inf && mpf_sgn(g->r.f)>0){
		mpf_set_si(a->r.f,1);
		a->r.inf = 1;
		a->r.open = g->r.open;
	}else if (g->r.inf && mpf_sgn(g->r.f)<0){
		mpf_set_ui(a->r.f,0);
		a->r.inf = 0;
		a->r.open = g->r.open;
	}else{
		mpf_exp(a->r.f, g->r.f,(PBITS+4)/3.322);
		a->r.inf = 0;
		a->r.open = g->r.open;
	}
	
	gnum_clear(&t1);
	gnum_clear(&t2);
	/* check gb */
	#if 0
	{
		printf("check gb: ");print_gb(a);printf("\n");
	}
	#endif
}

void powg(gbnd_s *a, const gbnd_s *x, const gbnd_s *y){
	int nan = 0;
	gbnd_s gb1,gb2,gb3,gb4,tmpgb;

	// NaN returns NaN 
	if (x->nan || y->nan) {
		a->nan = 1;
		AOP2("NaN",a,gb,x,gb,plusg,y,gb);
		return;
	}
	// x<0 ^ interval
	if (mpf_sgn(x->l.f)<0 && !exact(&y->l,&y->r)) {
		a->nan = 1;
		AOP2("NaN",a,gb,x,gb,plusg,y,gb);
		return;
	}
	
	a->nan = 0;
	gbnd_init(&gb1);
	gbnd_init(&gb2);
	gbnd_init(&gb3);
	gbnd_init(&gb4);
	gbnd_init(&tmpgb);
	
    // ------------------
	// x.left ^ y. left
	// ------------------
	power(&gb1.l,&gb1.r ,&x->l,&y->l,&nan);
	// ------------------
	// x.left ^ y.right
	// ------------------
	power(&gb2.l,&gb2.r,&x->l,&y->r,&nan);
	// ------------------
	// x.right ^ y.left
	// ------------------
	power(&gb3.l,&gb3.r,&x->r,&y->l,&nan);
	// ------------------
	// x.right ^ y.right
	// ------------------
	power(&gb4.l,&gb4.r,&x->r,&y->r,&nan);

	a->nan = nan;
	if (a->nan) {
		mpf_set_si(a->l.f, 0); mpf_set_si(a->r.f, 0);
		a->l.inf = a->r.inf = 0;
		a->l.open = a->r.open = 1;
		AOP2("NaN",a,gb,x,gb,plusg,y,gb);
	}else{
		combinegb(&tmpgb,&gb1,&gb2);
		combinegb(&tmpgb,&tmpgb,&gb3);
		combinegb(&tmpgb,&tmpgb,&gb4);

		a->l = tmpgb.l;
		a->r = tmpgb.r;
		/* See if 0 is in the range */
		if (mpf_sgn(a->l.f)>0){
			if (spanszerogQ(x)) {
				mpf_set_si(a->l.f, 0);
				a->l.inf = 0;
				a->l.open = 0;
			}
		}
	}
	/* check gb */
	#if 0
	{
		printf("check gb: ");print_gb(a);printf("\n");
	}
	#endif
}

/* use the corresponding method to caluclate power */
void power(gnum_s *a, gnum_s *b, const gnum_s *x, const gnum_s *y, int *nan){
	/* x=0, y=0 -> Nan
			y>0 -> 0
	        y<0 -> Nan
	*/ 
	if (!mpf_sgn(x->f) && !x->open){ 
		if (!mpf_sgn(y->f) && !y->open){*nan = 1;} 
		else if (mpf_sgn(y->f)>0){
			mpf_set_ui(a->f,0);
			a->inf = 0;
			a->open = 0;
		} else{*nan =1;} 
	/* x = +inf,  y=0 -> 1
			      y>0 -> +inf
	        	  y<0 -> 0
	*/ 
	}else if (mpf_sgn(x->f)>0 && x->inf){ 
		if (!mpf_sgn(y->f) && !y->open){
			mpf_set_ui(a->f,1);
			a->inf = 0;
			a->open = 0;
		} else if (mpf_sgn(y->f)>0){
			mpf_set_ui(a->f,1);
			a->inf = 1;
			a->open = x->open;
		} else{
			mpf_set_ui(a->f,0);
			a->inf = 0;
			a->open = x->open;
		} 
	/* x>0,  y=0 -> 1
	         y>0 -> +inf -> +inf
	 			 -> x ^ y
	         y<0 -> -inf -> 0
				 -> x ^ y
	*/ 
	}else if (mpf_sgn(x->f)>0){ 
		if (!mpf_sgn(y->f) && !y->open){
			mpf_set_ui(a->f,1);
			a->inf = 0;
			a->open = 0;
		} else if (mpf_sgn(y->f)>0){ 
			if (y->inf){
				mpf_set_ui(a->f,1);
				a->inf = 1;
				a->open = y->open;
			} else{ 
				mpf_pow(a->f,b->f,x->f,y->f,PBITS/3.322);
				a->inf = 0; 
				a->open = x->open || y->open;
				b->inf = 0; 
				b->open = x->open || y->open;
			}
		}else{ 
			if (y->inf){
				mpf_set_ui(a->f,0);
				a->inf = 0;
				a->open = y->open;
			}else{
				mpf_pow(a->f,b->f,x->f,y->f,(PBITS+4)/3.322);
				a->inf = 0; 
				a->open = x->open || y->open;
				b->inf = 0; 
				b->open = x->open || y->open;
			} 
		} 
	/* x=-inf,  y=0 -> 1
				y>0 -> inf -> Nan
	                -> x ^ y
	            y<0 -> -inf -> Nan
	                -> 0
	*/ 
	}else if(mpf_sgn(x->f)<0 && x->inf){ 
		if (!mpf_sgn(y->f) && !y->open){
			mpf_set_ui(a->f,1);
			a->inf = 0;
			a->open = 0;
		} else if (mpf_sgn(y->f)>0){ 
			if (y->inf){*nan = 1;} 
			else{ 
				mpf_pow(a->f,b->f,x->f,y->f,(PBITS+4)/3.322);
				a->inf = 1; 
				a->open = x->open;
				b->inf = 1; 
				b->open = x->open;
			} 
		}else{ 
			if (y->inf){*nan = 1;} 
			else{
				mpf_set_ui(a->f,0);
				a->inf = 0;
				a->open = a->open;} 
		} 
	/* x<0,  y=0 -> 1
		   	 y>0 -> inf -> Nan
	             -> x ^ y
	         y<0 -> -inf -> Nan
	             -> x ^ y
	*/ 
	}else{ 
		if (!mpf_sgn(y->f) && !y->open){
			mpf_set_ui(a->f,1);
			a->inf = 0;
			a->open = 0;
		} else if (mpf_sgn(y->f)>0){ 
			if (y->inf){*nan = 1;} 
			else{ 
				mpf_pow(a->f,b->f,x->f,y->f,(PBITS+4)/3.322);
				a->inf = 0; 
				a->open = x->open || y->open;
				b->inf = 0; 
				b->open = x->open || y->open;
			} 
		} else{ 
			if (y->inf){*nan = 1;} 
			else{ 
				mpf_pow(a->f,b->f,x->f,y->f,(PBITS+4)/3.322);
				a->inf = 0; 
				a->open = x->open || y->open;
				b->inf = 0; 
				b->open = x->open || y->open;
			} 
		} 
	}
}

/* combine two gbounds into a gbound which contains all the values */
void combinegb(gbnd_s *interval, const gbnd_s *x, const gbnd_s *y){
	/* Left endpoint of the interval */
	/* -Infinity case */
	if ((x->l.inf && mpf_sgn(x->l.f)<0) || (y->l.inf && mpf_sgn(y->l.f)<0)){
		if ((x->l.inf && mpf_sgn(x->l.f)<0) && (y->l.inf && mpf_sgn(y->l.f)<0)){
			mpf_set_si(interval->l.f,-1);
			interval->l.inf = 1;
			interval->l.open = x->l.open && y->l.open;
		}else{
			if (x->l.inf && mpf_sgn(x->l.f)<0){
				mpf_set_si(interval->l.f,-1);
				interval->l.inf = 1;
				interval->l.open = x->l.open;
			}else{
				mpf_set_si(interval->l.f,-1);
				interval->l.inf = 1;
				interval->l.open = y->l.open;
			}
		}
	/* Infinity case */
	}else if ((x->l.inf && mpf_sgn(x->l.f)>0) || (y->l.inf && mpf_sgn(y->l.f)>0)){
		if ((x->l.inf && mpf_sgn(x->l.f)>0) && (y->l.inf && mpf_sgn(y->l.f)>0)){
			mpf_set_si(interval->l.f,1);
			interval->l.inf = 1;
			interval->l.open = x->l.open || y->l.open;
		}else{
			if (x->l.inf && mpf_sgn(x->l.f)>0){
				mpf_set(interval->l.f,y->l.f);
				interval->l.inf = y->l.inf;
				interval->l.open = y->l.open;
			}else{
				mpf_set(interval->l.f,x->l.f);
				interval->l.inf = x->l.inf;
				interval->l.open = x->l.open;
			}
		}
	/* General case, find minimum */
	}else{
		if (mpf_cmp(x->l.f,y->l.f)>0){
			mpf_set(interval->l.f,y->l.f);
			interval->l.inf = y->l.inf;
			interval->l.open = y->l.open;
		}else if (mpf_cmp(x->l.f,y->l.f)<0){
			mpf_set(interval->l.f,x->l.f);
			interval->l.inf = x->l.inf;
			interval->l.open = x->l.open;
		}else{
			mpf_set(interval->l.f,x->l.f);
			interval->l.inf = x->l.inf;
			interval->l.open = x->l.open && y->l.open;
		}
	}

	/* Right endpoint of the interval */
	/* -Infinity case */
	if ((x->r.inf && mpf_sgn(x->r.f)<0) || (y->r.inf && mpf_sgn(y->r.f)<0)){
		if ((x->r.inf && mpf_sgn(x->r.f)<0) && (y->r.inf && mpf_sgn(y->r.f)<0)){
			mpf_set_si(interval->r.f,-1);
			interval->r.inf= 1;
			interval->r.open = x->r.open || y->r.open;
		}else{
			if (x->r.inf && mpf_sgn(x->r.f)<0){
				mpf_set(interval->r.f,y->r.f);
				interval->r.inf = y->r.inf;
				interval->r.open = y->r.open;
			}else{
				mpf_set(interval->r.f,x->r.f);
				interval->r.inf = x->r.inf;
				interval->r.open = x->r.open;
			}
		}
	/* Infinity case */
	}else if ((x->r.inf && mpf_sgn(x->r.f)>0) || (y->r.inf && mpf_sgn(y->r.f)>0)){
		if ((x->r.inf && mpf_sgn(x->r.f)>0) && (y->r.inf && mpf_sgn(y->r.f)>0)){
			mpf_set_si(interval->r.f,1);
			interval->r.inf = 1;
			interval->r.open = x->r.open && y->r.open;
		}else{
			if (x->r.inf && mpf_sgn(x->r.f)>0){
				mpf_set_ui(interval->r.f,1);
				interval->r.inf = 1;
				interval->r.open= x->r.open;
			}else{
				mpf_set_ui(interval->r.f,1);
				interval->r.inf = 1;
				interval->r.open = y->r.open;
			}
		}
	/* General case, find maximum */
	}else{
		if (mpf_cmp(x->r.f,y->r.f)<0){
			mpf_set(interval->r.f,y->r.f);
			interval->r.inf = y->r.inf;
			interval->r.open= y->r.open;
		}else if (mpf_cmp(x->r.f,y->r.f)>0){
			mpf_set(interval->r.f,x->r.f);
			interval->r.inf = x->r.inf;
			interval->r.open = x->r.open;
		}else{
			mpf_set(interval->r.f,x->r.f);
			interval->r.inf = x->r.inf;
			interval->r.open = x->r.open && y->r.open;
		}
	}
}

/* concat two char list */
void concat(char* result, const char *s1, const char *s2){
    strcpy(result, s1);
    strcat(result, s2);
}

// convert a.b to a/b
void dec2fra(fraction *a,const mpf_t b){
	mpf_t tmp;
	mpf_init(tmp);
	mpz_t num,den,gcd;
	mpz_init(num);
	mpz_init(den);
	mpz_init(gcd);

	mpf_abs(tmp,b);
	char top[256],bot[256];
	mp_exp_t exp;
	mpf_get_str (top, &exp, 10, 0, tmp);
	mpz_set_str(num,top,10);
	mpz_set_ui(den,10);
	mpz_pow_ui(den,den,strlen(top)-exp);
	mpz_abs(num,num);
	mpz_abs(den,den);
	mpz_gcd(gcd,num,den);
	mpz_div(num,num,gcd);
	mpz_div(den,den,gcd);
	mpz_set(a->num,num);
	mpz_set(a->den,den);
	a->positive = (mpf_sgn(b) == -1) ? 0 : 1;
	mpz_clear(num);
	mpz_clear(den);
	mpz_clear(gcd);
	mpf_clear(tmp);
}

/* return a string between begin and end */
void substring(char *res, char *str, size_t begin, size_t end) { 
	if (str == 0 || strlen(str) == 0 || strlen(str) < begin || strlen(str) < end || begin == end){
		res[0] = '0';
		res[1] = '\0';
	}else{
		strncpy(res,str+begin,end-begin);
		res[end-begin] = '\0';
	}
}

/* check whether gnum is a interval or exact number
	1 indicates exact, 0 indicates interval
*/
int exact(const gnum_s *a,const gnum_s *b){
	return !mpf_cmp(a->f,b->f) && a->inf==b->inf && a->open==b->open;
}

/* find the first extended precision term */
void smallest_term(mpf_t min, int prec){
	mp_bitcnt_t bits =  prec * 3.322 + 4;
	unsigned int smallest = prec * 3.322 + 1.0;
	mpf_t one;
	mpf_init2 (one,bits);
	mpf_set_ui (one, 1);
	/* right shifting 1 to the (last+1) bit in the precision */
	mpf_div_2exp (min, one, smallest);
	mpf_clear (one);
}

/* find Taylor Series for logarithm
	 x - (x^2)/2 + (x^3)/3 ...
 */
void Taylor_Series_Log(mpf_t log, mpf_t x, unsigned int prec){
	mp_bitcnt_t bits = prec * 3.322 + 4;
	mpf_t top, term;
	mpf_init2 (top, bits);
	mpf_init2 (term, bits);
	mpf_t minterm;
	mpf_init2 (minterm, bits+4);
	smallest_term (minterm, prec+1);

	mpf_set_ui(top, 1);
	mpf_set_ui(log, 0);
	int n=1;
	while (1){
		mpf_mul(top,top,x);
		mpf_set(term,top);
		mpf_div_ui (term, term, n);
		mpf_add (log, log, term);
		/* if the term is smaller than the minimum value, it does not 
		   affect the result, terminate from here
		*/
		mpf_abs (term, term);
		if (mpf_cmp (term, minterm)<0) break;

		n ++;
	}
	mpf_clear (term);
	mpf_clear (top);
	mpf_clear (minterm);
}

/* find Taylor Series for exponential
	 1 + x/1! + (x^2)/2! ...
*/
void Taylor_Series_Exp(mpf_t exp, const mpf_t z, unsigned int prec){
	mpf_t top, fact, term;

	mp_bitcnt_t bits = prec * 3.322+4;
	mpf_init2 (top, bits);
	mpf_init2 (fact, bits);
	mpf_init2 (term, bits);

	mpf_t minterm;
	mpf_init2 (minterm, bits+4);
	smallest_term (minterm, prec+1);

	mpf_set_ui (top, 1);
	mpf_set_ui (exp, 1);
	int n=1;
	while (1){
		mpf_mul(top,top,z);
		mpf_set(term,top);
		/* find n! */
		factorial(fact,n);
		mpf_div(term, term, fact);
		mpf_add (exp, exp, term);

		/* if the term is smaller than the minimum value, it does not 
		   affect the result, terminate from here
		*/
		mpf_abs (term, term);
		if (mpf_cmp (term, minterm)<0) break;

		n++;
	}
	mpf_clear (top);
	mpf_clear (fact);
	mpf_clear (term);
	mpf_clear (minterm);
}

/* calculate n! */
void factorial (mpf_t res, unsigned int n){
	mpz_t fact;
	/* find how many bits needed to calculate n!
		otherwise, if extend the bits n! = inf 
	*/
	unsigned int digits = find_fact_digits(n);
	mpz_init2 (fact,digits);
	mpz_fac_ui (fact, n);
	mpf_set_z (res, fact);

	mpz_clear (fact);
}

/* find how many digits in n!
	n! = sum from 1 to n (log_10(n)+1)
	   = (nln(n) - n - 1) / ln(10)
 */
unsigned int find_fact_digits(const int n){
	mpf_t tmp,one,ln10;
	mpf_init2(tmp,PBITS+4);
	mpf_init2(one,PBITS+4);
	mpf_init2(ln10,PBITS+4);
	mpf_set_str(ln10,"2.30258509299404568401799145468436420760110148862877297603332790096757260967735248023599720508959829834196778404228",10);

	/* (nln(n) - n - 1) / ln(10) */
	mpf_set_ui(tmp,n);
	mpf_log(tmp,tmp,(PBITS+4)/3.322);
	mpf_mul_ui(tmp,tmp,n);
	mpf_sub_ui(tmp,tmp,n);
	mpf_add_ui(tmp,tmp,1);
	mpf_div(tmp,tmp,ln10);
	/* add one extra bits as insurance and also remedy the loss in ln(10) */
	int digits = mpf_get_ui(tmp)+1;

	mpf_clear(one);
	mpf_clear(ln10);
	mpf_clear(tmp);
	return digits;
}

/* find how many digits in x^y
	x^y = log_10(x^y)
	    = y * ln(x) / ln(10)
*/
unsigned int find_pow_digits(const mpf_t x, const mpf_t y,unsigned int prec){
	mp_bitcnt_t bits = prec * 3.322;
	mpf_t tmp,ln10;
	mpz_t res;
	unsigned int result;
	mpf_init2(tmp,bits);
	mpf_init2(ln10,bits);
	mpz_init2(res,bits);
	mpf_set_str(ln10,"2.302585092994045684017991454684364207601101488628772976033",10);
	mpf_log(tmp,x,prec);
	mpf_mul(tmp,tmp,y);
	mpf_div(tmp,tmp,ln10);
	mpf_ceil(tmp,tmp);
	mpz_set_f(res,tmp);
	result = mpz_get_ui(res);
	mpf_clear(tmp);
	mpf_clear(ln10);
	mpz_clear(res);
	return result;
}

/* extract one ulp base on the precision */
void ulp(mpf_t small, int prec){
	mp_bitcnt_t bits =  prec * 3.322;
	char smallest[20];
	mpf_t tmp;
	mpf_init2 (tmp,bits);
	sprintf(smallest,"1e-%i", prec);
	mpf_set_str (small, smallest,10);
	mpf_clear (tmp);
}

/* Test if interval g is not nowhere (somewhere) equal to interval h. */

int nneqgQ(const gbnd_s *g, const gbnd_s *h)
{
	return (g->nan && h->nan) || (!(g->nan || h->nan) && !(ltgQ(g, h) || gtgQ(g, h)));
}