/*
 * @BEGIN LICENSE
 *
 * ivocalc by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#include "psi4/libmints/vector.h" 
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libfock/jk.h"
#include "ivo.h"

namespace psi{ namespace ivo {

ivocalc::ivocalc(SharedWavefunction ref_wfn, Options& options): Wavefunction(options) {
	outfile->Printf("\n ------ Improved Virtual Orbital SETX ------ \n");
	shallow_copy(ref_wfn);
	reference_wavefunction_ = ref_wfn;

	print = options_.get_int("PRINT");
	hole = options_.get_int("Hole");

	//1. Get necessary information
	std::shared_ptr<PSIO> psio(_default_psio_lib_);
	if (!ref_wfn)
		throw PSIEXCEPTION("SCF has not been run yet!");

	outfile->Printf("\n C alpha ----------------------------------");
	Ca_->print();
}

ivocalc::~ivocalc() {}

double ivocalc::compute_energy() {
	outfile->Printf("%d", nirrep_);
        outfile->Printf("%d", nmopi_[0]);
	outfile->Printf("%d", doccpi_[0]);
        outfile->Printf("\n ");
	SharedMatrix C_occ(new Matrix("C_occ", nirrep_, nmopi_, doccpi_));
        C_occ->print();
	//1. Form Coeffs with a hole
	for (int h = 0; h < nirrep_; ++h) {
		for (int u = 0; u < nmopi_[h]; ++u) {
			for (int i = 0; i < doccpi_[h]; ++i) {
				if (i == hole) {
					C_occ->set(h, u, i, 0.0);
				}
				else {
					C_occ->set(h, u, i, Ca_->get(h, u, i));
				}
			}
		}
	}

	//2. cumpute G n-1 uv
	std::shared_ptr<JK> JK_occ;
	JK_occ = JK::build_JK(reference_wavefunction_->basisset(), reference_wavefunction_->get_basisset("DF_BASIS_SCF"), options_);
	JK_occ->set_memory(Process::environment.get_memory() * 0.8);

	// JK_core->set_cutoff(options_.get_double("INTEGRAL_SCREENING"));
	JK_occ->set_cutoff(options_.get_double("INTEGRAL_SCREENING"));
	JK_occ->initialize();
	JK_occ->set_do_J(true);
	// JK_core->set_allow_desymmetrization(true);
	JK_occ->set_do_K(true);

	std::vector<std::shared_ptr<Matrix>>& Cl = JK_occ->C_left();
	std::vector<std::shared_ptr<Matrix>>& Cr = JK_occ->C_right();
	Cl.clear();
	Cr.clear();
	Cl.push_back(C_occ);
	Cr.push_back(C_occ);

	JK_occ->compute();

	//3. Form Fuv n-1 alpha and beta

	SharedMatrix F_nm1 = JK_occ->J()[0];
	SharedMatrix K_nm1 = JK_occ->K()[0];

	F_nm1->scale(2.0);
	F_nm1->subtract(K_nm1); //2J-K

	SharedMatrix T = SharedMatrix(reference_wavefunction_->matrix_factory()->create_matrix(PSIF_SO_T)); //initialize matrix of wfn's SO size
	SharedMatrix V_oe = SharedMatrix(reference_wavefunction_->matrix_factory()->create_matrix(PSIF_SO_V));

	MintsHelper mints(reference_wavefunction_); //Mintshelper read many ints (T, V, Prop, ...)
	T = mints.so_kinetic(); //read from mints
	V_oe = mints.so_potential();

	//Now build alpha/beta oei in AO basis
	SharedMatrix Ha = T->clone();
	//SharedMatrix Hb = T->clone();

	Ha->add(V_oe);
	//Hb->add(V_oe);

	// transform h_pq to h_uv
	Ha->transform(Ca_);
	F_nm1->transform(Ca_);

	F_nm1->add(Ha); //F_uv = H_uv + G_uv

	if(print == 2){
		outfile->Printf("\n ** F_n-1_uv checkpoint ----------------------------------");
		F_nm1->print();
	}

	//4. Solve Fvv n-1 Uvv = epi_v * Uvv
	//Get F n-1 vv
	Dimension zeropi = doccpi_;
	dvirpi_ = nmopi_ - doccpi_;
	for (int h = 0; h < nirrep_; ++h) {
		zeropi[h] = nmopi_[h] - nmopi_[h];
	}
	Slice oo(zeropi, doccpi_);
	Slice vv(doccpi_, nmopi_);

	SharedMatrix F_nm1_vv = F_nm1->get_block(vv, vv);
	
	if(print == 2) {
		outfile->Printf("\n **Start Checkpoint**");
		F_nm1_vv->print();
		outfile->Printf("\n **End Checkpoint**");
	}

	//diagonalize F_n-1_vv

	SharedMatrix Uvv(new Matrix("U_vv", nirrep_, dvirpi_, dvirpi_));
	SharedVector epiv(new Vector("E_v", nirrep_, dvirpi_));

	F_nm1_vv->diagonalize(Uvv, epiv);
	
	if(print == 2){
		Uvv->print();
		epiv->print();
	}

	//5. C ivocalc = C ori * (I direct+ Uvv)
	SharedMatrix Utran(new Matrix("Iocc direct+ Uvv", nirrep_, nmopi_, nmopi_));

	SharedMatrix Uoo(new Matrix("Iocc", nirrep_, doccpi_, doccpi_));
	for (int h = 0; h < nirrep_; ++h) {
		for (int i = 0; i < doccpi_[h]; ++i) {
			Uoo->set(h, i, i, 1.0);
		}
	}

	Utran->set_block(oo, oo, Uoo);
	Utran->set_block(vv, vv, Uvv);

	if(print == 2){
		Utran->print();
	}

	SharedMatrix Civocalc(new Matrix("C_ivocalc", nirrep_, nmopi_, nmopi_));
	Civocalc->accumulate_product(Ca_, Utran);
	Civocalc->print();
	Ca_->copy(Civocalc);
	Cb_ = Ca_;

	return 0.0;
}

}} // End namespaces

