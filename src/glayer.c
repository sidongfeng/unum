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

#include "glayer.h"
#include "uenv.h" /* PBITS */

void gnum_init(gnum_s *gn)
{
	mpf_init2(gn->f, PBITS);
}

void gnum_clear(gnum_s *gn)
{
	mpf_clear(gn->f);
}

void gbnd_init(gbnd_s *g)
{
	mpf_init2(g->l.f, PBITS);
	mpf_init2(g->r.f, PBITS);
}

void gbnd_clear(gbnd_s *g)
{
	mpf_clear(g->l.f);
	mpf_clear(g->r.f);
}

void nth_root_init(nth_root_par *nrp)
{
	mpf_init2(nrp->B, PBITS);
	mpf_init2(nrp->alpha, PBITS);
	mpf_init2(nrp->r, PBITS);
	mpf_init2(nrp->y, PBITS);
	mpf_init2(nrp->y_l, PBITS);
	nrp->beta = 0;
}

void nth_root_clear(nth_root_par *nrp)
{
	mpf_clear(nrp->B);
	mpf_clear(nrp->alpha);
	mpf_clear(nrp->r);
	mpf_clear(nrp->y);
	mpf_clear(nrp->y_l);
}

void fraction_init(fraction *fra)
{
	mpz_init2(fra->num,PBITS);
	mpz_init2(fra->den,PBITS);
	fra->positive = 0;
}

void fraction_clear(fraction *fra)
{
	mpz_clear(fra->num);
	mpz_clear(fra->den);
}
