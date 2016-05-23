/**
 *  @file    pbamparm.c
 *  @ingroup PBAMparm
 *  @author  Andrew Stevens
 *  @brief   Class PBAMparm methods
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2014 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "pbamparm.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */


VPUBLIC PBAMparm* PBAMparm_ctor(PBAMparm_CalcType type) {

    /* Set up the structure */
    PBAMparm *thee = VNULL;
    thee = (PBAMparm*)Vmem_malloc(VNULL, 1, sizeof(PBAMparm));
    VASSERT( thee != VNULL);
    VASSERT( PBAMparm_ctor2(thee, type) == VRC_SUCCESS );

    return thee;
}

VPUBLIC Vrc_Codes PBAMparm_ctor2(PBAMparm *thee, PBAMparm_CalcType type) {

    int i;

    if (thee == VNULL) return VRC_FAILURE;

    thee->parsed = 0; //TODO: Initialization list taken from geoflow
    thee->type = type;
    thee->vdw = 0;
    thee->etol = 1.0e-6;

    return VRC_SUCCESS;
}

VPUBLIC void PBAMparm_dtor(PBAMparm **thee) {
    if ((*thee) != VNULL) {
        PBAMparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(PBAMparm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void PBAMparm_dtor2(PBAMparm *thee) { ; }

VPUBLIC Vrc_Codes PBAMparm_check(PBAMparm *thee) {

    Vrc_Codes rc;

    rc = VRC_SUCCESS;

    Vnm_print(0, "PBAMparm_check:  checking PBAMparm object of type %d.\n",
              thee->type);

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "PBAMparm_check:  not filled!\n");
        return VRC_FAILURE;
    }


    /* Check type settings */
    //if ((thee->type != GFCT_MANUAL)&& (thee->type != GFCT_AUTO)&& (thee->type != GFCT_NONE)) {
    if(thee->type != GFCT_AUTO) {
         Vnm_print(2,"PBAMparm_check: type not set");
         rc = VRC_FAILURE;
    }

    return rc;
}

VPUBLIC void PBAMparm_copy(PBAMparm *thee, PBAMparm *parm) {
    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    thee->type = parm->type;
    thee->parsed = parm->parsed;

    thee->vdw = parm->vdw;
    thee->etol = parm->etol;
}

Vrc_Codes FUBAR(const char* name){
    Vnm_print(2, "parsePBAM:  ran out of tokens on %s!\n", name);
    return VRC_WARNING;
}

Vrc_Codes parseNonNeg(double* tf, double def, int* set, char* name, Vio* sock){
    char tok[VMAX_BUFSIZE];
    if(Vio_scanf(sock, "%s", tok) == 0) {
        *tf = def;
        return FUBAR(name);
    }

    if (sscanf(tok, "%lf", tf) == 0){
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing %s keyword!\n", tok, name);
        *tf = def;
        return VRC_WARNING;
    }else if(*tf < 0.0){
        Vnm_print(2, "parsePBAM:  %s must be greater than 0!\n", name);
        *tf = def;
        return VRC_WARNING;
    }

    *set = 1;
    return VRC_SUCCESS;
}

VPRIVATE Vrc_Codes PBAMparm_parseVDW(PBAMparm *thee, Vio *sock){
    const char* name = "vdw";
    char tok[VMAX_BUFSIZE];
	int tf;
    if(Vio_scanf(sock, "%s", tok) == 0) {
        return FUBAR(name);
    }


    if (sscanf(tok, "%u", &tf) == 0){
        Vnm_print(2, "NOsh:  Read non-unsigned int (%s) while parsing %s keyword!\n", tok, name);
        return VRC_WARNING;
    }else if(tf != 0 && tf != 1){
        Vnm_print(2, "parsePBAM:  %s must be 0 or 1!\n", name);
        return VRC_WARNING;
    }else{
        thee->vdw = tf;
    }
    thee->setvdw = 1;
    return VRC_SUCCESS;
}

VPRIVATE Vrc_Codes PBAMparm_parseETOL(PBAMparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing etol keyword!\n", tok);
		return VRC_WARNING;
	} else if(tf <= 0.0) {
		Vnm_print(2,"parsePBAM: etol must be greater than 0!\n");
		return VRC_WARNING;
	} else {
		thee->etol = tf;
	}


	return VRC_SUCCESS;


	VERROR1:
		Vnm_print(2, "parsePBAM: ran out of tokens!\n");
		return VRC_WARNING;

}

VPUBLIC Vrc_Codes PBAMparm_parseToken(PBAMparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock) {

    if (thee == VNULL) {
        Vnm_print(2, "parsePBAM:  got NULL thee!\n");
        return VRC_WARNING;
    }
    if (sock == VNULL) {
        Vnm_print(2, "parsePBAM:  got NULL socket!\n");
        return VRC_WARNING;
    }

    Vnm_print(0, "PBAMparm_parseToken:  trying %s...\n", tok);

    if (Vstring_strcasecmp(tok, "vdwdisp") == 0) {
        return PBAMparm_parseVDW(thee, sock);
    } else if(Vstring_strcasecmp(tok, "etol") == 0){
    	return PBAMparm_parseETOL(thee, sock);
    }else {
        Vnm_print(2, "parsePBAM:  Unrecognized keyword (%s)!\n", tok);
        return VRC_WARNING;
    }


    return VRC_FAILURE;
}
