void flexmat6::update_as_pqr_stu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vr*iNp*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vu*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqr_stu(){
    if(Npqr_stu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vr*iNq*iNp;
        locations.col(1) = vs + vt*iNs + vu*iNt*iNs;
        Vpqr_stu = sp_mat(locations.t(), vValues, iNp*iNq*iNr,iNs*iNt*iNu);
        Npqr_stu = 1;
        return Vpqr_stu;
    }
    else{
        return Vpqr_stu;
    }
}

void flexmat6::update_as_pqr_sut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vr*iNp*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vt*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqr_sut(){
    if(Npqr_sut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vr*iNq*iNp;
        locations.col(1) = vs + vu*iNs + vt*iNu*iNs;
        Vpqr_sut = sp_mat(locations.t(), vValues, iNp*iNq*iNr,iNs*iNu*iNt);
        Npqr_sut = 1;
        return Vpqr_sut;
    }
    else{
        return Vpqr_sut;
    }
}

void flexmat6::update_as_pqr_tsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vr*iNp*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vu*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqr_tsu(){
    if(Npqr_tsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vr*iNq*iNp;
        locations.col(1) = vt + vs*iNt + vu*iNs*iNt;
        Vpqr_tsu = sp_mat(locations.t(), vValues, iNp*iNq*iNr,iNt*iNs*iNu);
        Npqr_tsu = 1;
        return Vpqr_tsu;
    }
    else{
        return Vpqr_tsu;
    }
}

void flexmat6::update_as_pqr_tus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vr*iNp*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vs*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqr_tus(){
    if(Npqr_tus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vr*iNq*iNp;
        locations.col(1) = vt + vu*iNt + vs*iNu*iNt;
        Vpqr_tus = sp_mat(locations.t(), vValues, iNp*iNq*iNr,iNt*iNu*iNs);
        Npqr_tus = 1;
        return Vpqr_tus;
    }
    else{
        return Vpqr_tus;
    }
}

void flexmat6::update_as_pqr_ust(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vr*iNp*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vt*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqr_ust(){
    if(Npqr_ust == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vr*iNq*iNp;
        locations.col(1) = vu + vs*iNu + vt*iNs*iNu;
        Vpqr_ust = sp_mat(locations.t(), vValues, iNp*iNq*iNr,iNu*iNs*iNt);
        Npqr_ust = 1;
        return Vpqr_ust;
    }
    else{
        return Vpqr_ust;
    }
}

void flexmat6::update_as_pqr_uts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vr*iNp*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vs*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqr_uts(){
    if(Npqr_uts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vr*iNq*iNp;
        locations.col(1) = vu + vt*iNu + vs*iNt*iNu;
        Vpqr_uts = sp_mat(locations.t(), vValues, iNp*iNq*iNr,iNu*iNt*iNs);
        Npqr_uts = 1;
        return Vpqr_uts;
    }
    else{
        return Vpqr_uts;
    }
}

void flexmat6::update_as_pqs_rtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vs*iNp*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vu*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqs_rtu(){
    if(Npqs_rtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vs*iNq*iNp;
        locations.col(1) = vr + vt*iNr + vu*iNt*iNr;
        Vpqs_rtu = sp_mat(locations.t(), vValues, iNp*iNq*iNs,iNr*iNt*iNu);
        Npqs_rtu = 1;
        return Vpqs_rtu;
    }
    else{
        return Vpqs_rtu;
    }
}

void flexmat6::update_as_pqs_rut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vs*iNp*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vt*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqs_rut(){
    if(Npqs_rut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vs*iNq*iNp;
        locations.col(1) = vr + vu*iNr + vt*iNu*iNr;
        Vpqs_rut = sp_mat(locations.t(), vValues, iNp*iNq*iNs,iNr*iNu*iNt);
        Npqs_rut = 1;
        return Vpqs_rut;
    }
    else{
        return Vpqs_rut;
    }
}

void flexmat6::update_as_pqs_tru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vs*iNp*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vu*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqs_tru(){
    if(Npqs_tru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vs*iNq*iNp;
        locations.col(1) = vt + vr*iNt + vu*iNr*iNt;
        Vpqs_tru = sp_mat(locations.t(), vValues, iNp*iNq*iNs,iNt*iNr*iNu);
        Npqs_tru = 1;
        return Vpqs_tru;
    }
    else{
        return Vpqs_tru;
    }
}

void flexmat6::update_as_pqs_tur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vs*iNp*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vr*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqs_tur(){
    if(Npqs_tur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vs*iNq*iNp;
        locations.col(1) = vt + vu*iNt + vr*iNu*iNt;
        Vpqs_tur = sp_mat(locations.t(), vValues, iNp*iNq*iNs,iNt*iNu*iNr);
        Npqs_tur = 1;
        return Vpqs_tur;
    }
    else{
        return Vpqs_tur;
    }
}

void flexmat6::update_as_pqs_urt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vs*iNp*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vt*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqs_urt(){
    if(Npqs_urt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vs*iNq*iNp;
        locations.col(1) = vu + vr*iNu + vt*iNr*iNu;
        Vpqs_urt = sp_mat(locations.t(), vValues, iNp*iNq*iNs,iNu*iNr*iNt);
        Npqs_urt = 1;
        return Vpqs_urt;
    }
    else{
        return Vpqs_urt;
    }
}

void flexmat6::update_as_pqs_utr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vs*iNp*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vr*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqs_utr(){
    if(Npqs_utr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vs*iNq*iNp;
        locations.col(1) = vu + vt*iNu + vr*iNt*iNu;
        Vpqs_utr = sp_mat(locations.t(), vValues, iNp*iNq*iNs,iNu*iNt*iNr);
        Npqs_utr = 1;
        return Vpqs_utr;
    }
    else{
        return Vpqs_utr;
    }
}

void flexmat6::update_as_pqt_rsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vt*iNp*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vu*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqt_rsu(){
    if(Npqt_rsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vt*iNq*iNp;
        locations.col(1) = vr + vs*iNr + vu*iNs*iNr;
        Vpqt_rsu = sp_mat(locations.t(), vValues, iNp*iNq*iNt,iNr*iNs*iNu);
        Npqt_rsu = 1;
        return Vpqt_rsu;
    }
    else{
        return Vpqt_rsu;
    }
}

void flexmat6::update_as_pqt_rus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vt*iNp*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vs*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqt_rus(){
    if(Npqt_rus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vt*iNq*iNp;
        locations.col(1) = vr + vu*iNr + vs*iNu*iNr;
        Vpqt_rus = sp_mat(locations.t(), vValues, iNp*iNq*iNt,iNr*iNu*iNs);
        Npqt_rus = 1;
        return Vpqt_rus;
    }
    else{
        return Vpqt_rus;
    }
}

void flexmat6::update_as_pqt_sru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vt*iNp*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vu*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqt_sru(){
    if(Npqt_sru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vt*iNq*iNp;
        locations.col(1) = vs + vr*iNs + vu*iNr*iNs;
        Vpqt_sru = sp_mat(locations.t(), vValues, iNp*iNq*iNt,iNs*iNr*iNu);
        Npqt_sru = 1;
        return Vpqt_sru;
    }
    else{
        return Vpqt_sru;
    }
}

void flexmat6::update_as_pqt_sur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vt*iNp*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vr*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqt_sur(){
    if(Npqt_sur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vt*iNq*iNp;
        locations.col(1) = vs + vu*iNs + vr*iNu*iNs;
        Vpqt_sur = sp_mat(locations.t(), vValues, iNp*iNq*iNt,iNs*iNu*iNr);
        Npqt_sur = 1;
        return Vpqt_sur;
    }
    else{
        return Vpqt_sur;
    }
}

void flexmat6::update_as_pqt_urs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vt*iNp*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vs*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqt_urs(){
    if(Npqt_urs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vt*iNq*iNp;
        locations.col(1) = vu + vr*iNu + vs*iNr*iNu;
        Vpqt_urs = sp_mat(locations.t(), vValues, iNp*iNq*iNt,iNu*iNr*iNs);
        Npqt_urs = 1;
        return Vpqt_urs;
    }
    else{
        return Vpqt_urs;
    }
}

void flexmat6::update_as_pqt_usr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vt*iNp*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vr*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqt_usr(){
    if(Npqt_usr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vt*iNq*iNp;
        locations.col(1) = vu + vs*iNu + vr*iNs*iNu;
        Vpqt_usr = sp_mat(locations.t(), vValues, iNp*iNq*iNt,iNu*iNs*iNr);
        Npqt_usr = 1;
        return Vpqt_usr;
    }
    else{
        return Vpqt_usr;
    }
}

void flexmat6::update_as_pqu_rst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vu*iNp*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vt*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqu_rst(){
    if(Npqu_rst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vu*iNq*iNp;
        locations.col(1) = vr + vs*iNr + vt*iNs*iNr;
        Vpqu_rst = sp_mat(locations.t(), vValues, iNp*iNq*iNu,iNr*iNs*iNt);
        Npqu_rst = 1;
        return Vpqu_rst;
    }
    else{
        return Vpqu_rst;
    }
}

void flexmat6::update_as_pqu_rts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vu*iNp*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vs*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqu_rts(){
    if(Npqu_rts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vu*iNq*iNp;
        locations.col(1) = vr + vt*iNr + vs*iNt*iNr;
        Vpqu_rts = sp_mat(locations.t(), vValues, iNp*iNq*iNu,iNr*iNt*iNs);
        Npqu_rts = 1;
        return Vpqu_rts;
    }
    else{
        return Vpqu_rts;
    }
}

void flexmat6::update_as_pqu_srt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vu*iNp*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vt*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqu_srt(){
    if(Npqu_srt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vu*iNq*iNp;
        locations.col(1) = vs + vr*iNs + vt*iNr*iNs;
        Vpqu_srt = sp_mat(locations.t(), vValues, iNp*iNq*iNu,iNs*iNr*iNt);
        Npqu_srt = 1;
        return Vpqu_srt;
    }
    else{
        return Vpqu_srt;
    }
}

void flexmat6::update_as_pqu_str(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vu*iNp*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vr*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqu_str(){
    if(Npqu_str == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vu*iNq*iNp;
        locations.col(1) = vs + vt*iNs + vr*iNt*iNs;
        Vpqu_str = sp_mat(locations.t(), vValues, iNp*iNq*iNu,iNs*iNt*iNr);
        Npqu_str = 1;
        return Vpqu_str;
    }
    else{
        return Vpqu_str;
    }
}

void flexmat6::update_as_pqu_trs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vu*iNp*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vs*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqu_trs(){
    if(Npqu_trs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vu*iNq*iNp;
        locations.col(1) = vt + vr*iNt + vs*iNr*iNt;
        Vpqu_trs = sp_mat(locations.t(), vValues, iNp*iNq*iNu,iNt*iNr*iNs);
        Npqu_trs = 1;
        return Vpqu_trs;
    }
    else{
        return Vpqu_trs;
    }
}

void flexmat6::update_as_pqu_tsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNp - vu*iNp*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vr*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pqu_tsr(){
    if(Npqu_tsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vq*iNp + vu*iNq*iNp;
        locations.col(1) = vt + vs*iNt + vr*iNs*iNt;
        Vpqu_tsr = sp_mat(locations.t(), vValues, iNp*iNq*iNu,iNt*iNs*iNr);
        Npqu_tsr = 1;
        return Vpqu_tsr;
    }
    else{
        return Vpqu_tsr;
    }
}

void flexmat6::update_as_prq_stu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vq*iNp*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vu*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prq_stu(){
    if(Nprq_stu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vq*iNr*iNp;
        locations.col(1) = vs + vt*iNs + vu*iNt*iNs;
        Vprq_stu = sp_mat(locations.t(), vValues, iNp*iNr*iNq,iNs*iNt*iNu);
        Nprq_stu = 1;
        return Vprq_stu;
    }
    else{
        return Vprq_stu;
    }
}

void flexmat6::update_as_prq_sut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vq*iNp*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vt*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prq_sut(){
    if(Nprq_sut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vq*iNr*iNp;
        locations.col(1) = vs + vu*iNs + vt*iNu*iNs;
        Vprq_sut = sp_mat(locations.t(), vValues, iNp*iNr*iNq,iNs*iNu*iNt);
        Nprq_sut = 1;
        return Vprq_sut;
    }
    else{
        return Vprq_sut;
    }
}

void flexmat6::update_as_prq_tsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vq*iNp*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vu*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prq_tsu(){
    if(Nprq_tsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vq*iNr*iNp;
        locations.col(1) = vt + vs*iNt + vu*iNs*iNt;
        Vprq_tsu = sp_mat(locations.t(), vValues, iNp*iNr*iNq,iNt*iNs*iNu);
        Nprq_tsu = 1;
        return Vprq_tsu;
    }
    else{
        return Vprq_tsu;
    }
}

void flexmat6::update_as_prq_tus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vq*iNp*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vs*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prq_tus(){
    if(Nprq_tus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vq*iNr*iNp;
        locations.col(1) = vt + vu*iNt + vs*iNu*iNt;
        Vprq_tus = sp_mat(locations.t(), vValues, iNp*iNr*iNq,iNt*iNu*iNs);
        Nprq_tus = 1;
        return Vprq_tus;
    }
    else{
        return Vprq_tus;
    }
}

void flexmat6::update_as_prq_ust(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vq*iNp*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vt*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prq_ust(){
    if(Nprq_ust == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vq*iNr*iNp;
        locations.col(1) = vu + vs*iNu + vt*iNs*iNu;
        Vprq_ust = sp_mat(locations.t(), vValues, iNp*iNr*iNq,iNu*iNs*iNt);
        Nprq_ust = 1;
        return Vprq_ust;
    }
    else{
        return Vprq_ust;
    }
}

void flexmat6::update_as_prq_uts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vq*iNp*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vs*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prq_uts(){
    if(Nprq_uts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vq*iNr*iNp;
        locations.col(1) = vu + vt*iNu + vs*iNt*iNu;
        Vprq_uts = sp_mat(locations.t(), vValues, iNp*iNr*iNq,iNu*iNt*iNs);
        Nprq_uts = 1;
        return Vprq_uts;
    }
    else{
        return Vprq_uts;
    }
}

void flexmat6::update_as_prs_qtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vs*iNp*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vu*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prs_qtu(){
    if(Nprs_qtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vs*iNr*iNp;
        locations.col(1) = vq + vt*iNq + vu*iNt*iNq;
        Vprs_qtu = sp_mat(locations.t(), vValues, iNp*iNr*iNs,iNq*iNt*iNu);
        Nprs_qtu = 1;
        return Vprs_qtu;
    }
    else{
        return Vprs_qtu;
    }
}

void flexmat6::update_as_prs_qut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vs*iNp*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vt*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prs_qut(){
    if(Nprs_qut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vs*iNr*iNp;
        locations.col(1) = vq + vu*iNq + vt*iNu*iNq;
        Vprs_qut = sp_mat(locations.t(), vValues, iNp*iNr*iNs,iNq*iNu*iNt);
        Nprs_qut = 1;
        return Vprs_qut;
    }
    else{
        return Vprs_qut;
    }
}

void flexmat6::update_as_prs_tqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vs*iNp*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vu*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prs_tqu(){
    if(Nprs_tqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vs*iNr*iNp;
        locations.col(1) = vt + vq*iNt + vu*iNq*iNt;
        Vprs_tqu = sp_mat(locations.t(), vValues, iNp*iNr*iNs,iNt*iNq*iNu);
        Nprs_tqu = 1;
        return Vprs_tqu;
    }
    else{
        return Vprs_tqu;
    }
}

void flexmat6::update_as_prs_tuq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vs*iNp*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vq*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prs_tuq(){
    if(Nprs_tuq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vs*iNr*iNp;
        locations.col(1) = vt + vu*iNt + vq*iNu*iNt;
        Vprs_tuq = sp_mat(locations.t(), vValues, iNp*iNr*iNs,iNt*iNu*iNq);
        Nprs_tuq = 1;
        return Vprs_tuq;
    }
    else{
        return Vprs_tuq;
    }
}

void flexmat6::update_as_prs_uqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vs*iNp*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vt*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prs_uqt(){
    if(Nprs_uqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vs*iNr*iNp;
        locations.col(1) = vu + vq*iNu + vt*iNq*iNu;
        Vprs_uqt = sp_mat(locations.t(), vValues, iNp*iNr*iNs,iNu*iNq*iNt);
        Nprs_uqt = 1;
        return Vprs_uqt;
    }
    else{
        return Vprs_uqt;
    }
}

void flexmat6::update_as_prs_utq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vs*iNp*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vq*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prs_utq(){
    if(Nprs_utq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vs*iNr*iNp;
        locations.col(1) = vu + vt*iNu + vq*iNt*iNu;
        Vprs_utq = sp_mat(locations.t(), vValues, iNp*iNr*iNs,iNu*iNt*iNq);
        Nprs_utq = 1;
        return Vprs_utq;
    }
    else{
        return Vprs_utq;
    }
}

void flexmat6::update_as_prt_qsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vt*iNp*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vu*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prt_qsu(){
    if(Nprt_qsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vt*iNr*iNp;
        locations.col(1) = vq + vs*iNq + vu*iNs*iNq;
        Vprt_qsu = sp_mat(locations.t(), vValues, iNp*iNr*iNt,iNq*iNs*iNu);
        Nprt_qsu = 1;
        return Vprt_qsu;
    }
    else{
        return Vprt_qsu;
    }
}

void flexmat6::update_as_prt_qus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vt*iNp*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vs*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prt_qus(){
    if(Nprt_qus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vt*iNr*iNp;
        locations.col(1) = vq + vu*iNq + vs*iNu*iNq;
        Vprt_qus = sp_mat(locations.t(), vValues, iNp*iNr*iNt,iNq*iNu*iNs);
        Nprt_qus = 1;
        return Vprt_qus;
    }
    else{
        return Vprt_qus;
    }
}

void flexmat6::update_as_prt_squ(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vt*iNp*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vu*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prt_squ(){
    if(Nprt_squ == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vt*iNr*iNp;
        locations.col(1) = vs + vq*iNs + vu*iNq*iNs;
        Vprt_squ = sp_mat(locations.t(), vValues, iNp*iNr*iNt,iNs*iNq*iNu);
        Nprt_squ = 1;
        return Vprt_squ;
    }
    else{
        return Vprt_squ;
    }
}

void flexmat6::update_as_prt_suq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vt*iNp*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vq*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prt_suq(){
    if(Nprt_suq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vt*iNr*iNp;
        locations.col(1) = vs + vu*iNs + vq*iNu*iNs;
        Vprt_suq = sp_mat(locations.t(), vValues, iNp*iNr*iNt,iNs*iNu*iNq);
        Nprt_suq = 1;
        return Vprt_suq;
    }
    else{
        return Vprt_suq;
    }
}

void flexmat6::update_as_prt_uqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vt*iNp*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vs*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prt_uqs(){
    if(Nprt_uqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vt*iNr*iNp;
        locations.col(1) = vu + vq*iNu + vs*iNq*iNu;
        Vprt_uqs = sp_mat(locations.t(), vValues, iNp*iNr*iNt,iNu*iNq*iNs);
        Nprt_uqs = 1;
        return Vprt_uqs;
    }
    else{
        return Vprt_uqs;
    }
}

void flexmat6::update_as_prt_usq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vt*iNp*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vq*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::prt_usq(){
    if(Nprt_usq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vt*iNr*iNp;
        locations.col(1) = vu + vs*iNu + vq*iNs*iNu;
        Vprt_usq = sp_mat(locations.t(), vValues, iNp*iNr*iNt,iNu*iNs*iNq);
        Nprt_usq = 1;
        return Vprt_usq;
    }
    else{
        return Vprt_usq;
    }
}

void flexmat6::update_as_pru_qst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vu*iNp*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vt*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pru_qst(){
    if(Npru_qst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vu*iNr*iNp;
        locations.col(1) = vq + vs*iNq + vt*iNs*iNq;
        Vpru_qst = sp_mat(locations.t(), vValues, iNp*iNr*iNu,iNq*iNs*iNt);
        Npru_qst = 1;
        return Vpru_qst;
    }
    else{
        return Vpru_qst;
    }
}

void flexmat6::update_as_pru_qts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vu*iNp*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vs*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pru_qts(){
    if(Npru_qts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vu*iNr*iNp;
        locations.col(1) = vq + vt*iNq + vs*iNt*iNq;
        Vpru_qts = sp_mat(locations.t(), vValues, iNp*iNr*iNu,iNq*iNt*iNs);
        Npru_qts = 1;
        return Vpru_qts;
    }
    else{
        return Vpru_qts;
    }
}

void flexmat6::update_as_pru_sqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vu*iNp*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vt*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pru_sqt(){
    if(Npru_sqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vu*iNr*iNp;
        locations.col(1) = vs + vq*iNs + vt*iNq*iNs;
        Vpru_sqt = sp_mat(locations.t(), vValues, iNp*iNr*iNu,iNs*iNq*iNt);
        Npru_sqt = 1;
        return Vpru_sqt;
    }
    else{
        return Vpru_sqt;
    }
}

void flexmat6::update_as_pru_stq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vu*iNp*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vq*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pru_stq(){
    if(Npru_stq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vu*iNr*iNp;
        locations.col(1) = vs + vt*iNs + vq*iNt*iNs;
        Vpru_stq = sp_mat(locations.t(), vValues, iNp*iNr*iNu,iNs*iNt*iNq);
        Npru_stq = 1;
        return Vpru_stq;
    }
    else{
        return Vpru_stq;
    }
}

void flexmat6::update_as_pru_tqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vu*iNp*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vs*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pru_tqs(){
    if(Npru_tqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vu*iNr*iNp;
        locations.col(1) = vt + vq*iNt + vs*iNq*iNt;
        Vpru_tqs = sp_mat(locations.t(), vValues, iNp*iNr*iNu,iNt*iNq*iNs);
        Npru_tqs = 1;
        return Vpru_tqs;
    }
    else{
        return Vpru_tqs;
    }
}

void flexmat6::update_as_pru_tsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNp - vu*iNp*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vq*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pru_tsq(){
    if(Npru_tsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vr*iNp + vu*iNr*iNp;
        locations.col(1) = vt + vs*iNt + vq*iNs*iNt;
        Vpru_tsq = sp_mat(locations.t(), vValues, iNp*iNr*iNu,iNt*iNs*iNq);
        Npru_tsq = 1;
        return Vpru_tsq;
    }
    else{
        return Vpru_tsq;
    }
}

void flexmat6::update_as_psq_rtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vq*iNp*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vu*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psq_rtu(){
    if(Npsq_rtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vq*iNs*iNp;
        locations.col(1) = vr + vt*iNr + vu*iNt*iNr;
        Vpsq_rtu = sp_mat(locations.t(), vValues, iNp*iNs*iNq,iNr*iNt*iNu);
        Npsq_rtu = 1;
        return Vpsq_rtu;
    }
    else{
        return Vpsq_rtu;
    }
}

void flexmat6::update_as_psq_rut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vq*iNp*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vt*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psq_rut(){
    if(Npsq_rut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vq*iNs*iNp;
        locations.col(1) = vr + vu*iNr + vt*iNu*iNr;
        Vpsq_rut = sp_mat(locations.t(), vValues, iNp*iNs*iNq,iNr*iNu*iNt);
        Npsq_rut = 1;
        return Vpsq_rut;
    }
    else{
        return Vpsq_rut;
    }
}

void flexmat6::update_as_psq_tru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vq*iNp*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vu*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psq_tru(){
    if(Npsq_tru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vq*iNs*iNp;
        locations.col(1) = vt + vr*iNt + vu*iNr*iNt;
        Vpsq_tru = sp_mat(locations.t(), vValues, iNp*iNs*iNq,iNt*iNr*iNu);
        Npsq_tru = 1;
        return Vpsq_tru;
    }
    else{
        return Vpsq_tru;
    }
}

void flexmat6::update_as_psq_tur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vq*iNp*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vr*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psq_tur(){
    if(Npsq_tur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vq*iNs*iNp;
        locations.col(1) = vt + vu*iNt + vr*iNu*iNt;
        Vpsq_tur = sp_mat(locations.t(), vValues, iNp*iNs*iNq,iNt*iNu*iNr);
        Npsq_tur = 1;
        return Vpsq_tur;
    }
    else{
        return Vpsq_tur;
    }
}

void flexmat6::update_as_psq_urt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vq*iNp*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vt*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psq_urt(){
    if(Npsq_urt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vq*iNs*iNp;
        locations.col(1) = vu + vr*iNu + vt*iNr*iNu;
        Vpsq_urt = sp_mat(locations.t(), vValues, iNp*iNs*iNq,iNu*iNr*iNt);
        Npsq_urt = 1;
        return Vpsq_urt;
    }
    else{
        return Vpsq_urt;
    }
}

void flexmat6::update_as_psq_utr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vq*iNp*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vr*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psq_utr(){
    if(Npsq_utr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vq*iNs*iNp;
        locations.col(1) = vu + vt*iNu + vr*iNt*iNu;
        Vpsq_utr = sp_mat(locations.t(), vValues, iNp*iNs*iNq,iNu*iNt*iNr);
        Npsq_utr = 1;
        return Vpsq_utr;
    }
    else{
        return Vpsq_utr;
    }
}

void flexmat6::update_as_psr_qtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vr*iNp*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vu*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psr_qtu(){
    if(Npsr_qtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vr*iNs*iNp;
        locations.col(1) = vq + vt*iNq + vu*iNt*iNq;
        Vpsr_qtu = sp_mat(locations.t(), vValues, iNp*iNs*iNr,iNq*iNt*iNu);
        Npsr_qtu = 1;
        return Vpsr_qtu;
    }
    else{
        return Vpsr_qtu;
    }
}

void flexmat6::update_as_psr_qut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vr*iNp*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vt*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psr_qut(){
    if(Npsr_qut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vr*iNs*iNp;
        locations.col(1) = vq + vu*iNq + vt*iNu*iNq;
        Vpsr_qut = sp_mat(locations.t(), vValues, iNp*iNs*iNr,iNq*iNu*iNt);
        Npsr_qut = 1;
        return Vpsr_qut;
    }
    else{
        return Vpsr_qut;
    }
}

void flexmat6::update_as_psr_tqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vr*iNp*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vu*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psr_tqu(){
    if(Npsr_tqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vr*iNs*iNp;
        locations.col(1) = vt + vq*iNt + vu*iNq*iNt;
        Vpsr_tqu = sp_mat(locations.t(), vValues, iNp*iNs*iNr,iNt*iNq*iNu);
        Npsr_tqu = 1;
        return Vpsr_tqu;
    }
    else{
        return Vpsr_tqu;
    }
}

void flexmat6::update_as_psr_tuq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vr*iNp*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vq*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psr_tuq(){
    if(Npsr_tuq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vr*iNs*iNp;
        locations.col(1) = vt + vu*iNt + vq*iNu*iNt;
        Vpsr_tuq = sp_mat(locations.t(), vValues, iNp*iNs*iNr,iNt*iNu*iNq);
        Npsr_tuq = 1;
        return Vpsr_tuq;
    }
    else{
        return Vpsr_tuq;
    }
}

void flexmat6::update_as_psr_uqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vr*iNp*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vt*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psr_uqt(){
    if(Npsr_uqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vr*iNs*iNp;
        locations.col(1) = vu + vq*iNu + vt*iNq*iNu;
        Vpsr_uqt = sp_mat(locations.t(), vValues, iNp*iNs*iNr,iNu*iNq*iNt);
        Npsr_uqt = 1;
        return Vpsr_uqt;
    }
    else{
        return Vpsr_uqt;
    }
}

void flexmat6::update_as_psr_utq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vr*iNp*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vq*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psr_utq(){
    if(Npsr_utq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vr*iNs*iNp;
        locations.col(1) = vu + vt*iNu + vq*iNt*iNu;
        Vpsr_utq = sp_mat(locations.t(), vValues, iNp*iNs*iNr,iNu*iNt*iNq);
        Npsr_utq = 1;
        return Vpsr_utq;
    }
    else{
        return Vpsr_utq;
    }
}

void flexmat6::update_as_pst_qru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vt*iNp*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vu*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pst_qru(){
    if(Npst_qru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vt*iNs*iNp;
        locations.col(1) = vq + vr*iNq + vu*iNr*iNq;
        Vpst_qru = sp_mat(locations.t(), vValues, iNp*iNs*iNt,iNq*iNr*iNu);
        Npst_qru = 1;
        return Vpst_qru;
    }
    else{
        return Vpst_qru;
    }
}

void flexmat6::update_as_pst_qur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vt*iNp*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vr*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pst_qur(){
    if(Npst_qur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vt*iNs*iNp;
        locations.col(1) = vq + vu*iNq + vr*iNu*iNq;
        Vpst_qur = sp_mat(locations.t(), vValues, iNp*iNs*iNt,iNq*iNu*iNr);
        Npst_qur = 1;
        return Vpst_qur;
    }
    else{
        return Vpst_qur;
    }
}

void flexmat6::update_as_pst_rqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vt*iNp*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vu*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pst_rqu(){
    if(Npst_rqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vt*iNs*iNp;
        locations.col(1) = vr + vq*iNr + vu*iNq*iNr;
        Vpst_rqu = sp_mat(locations.t(), vValues, iNp*iNs*iNt,iNr*iNq*iNu);
        Npst_rqu = 1;
        return Vpst_rqu;
    }
    else{
        return Vpst_rqu;
    }
}

void flexmat6::update_as_pst_ruq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vt*iNp*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vq*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pst_ruq(){
    if(Npst_ruq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vt*iNs*iNp;
        locations.col(1) = vr + vu*iNr + vq*iNu*iNr;
        Vpst_ruq = sp_mat(locations.t(), vValues, iNp*iNs*iNt,iNr*iNu*iNq);
        Npst_ruq = 1;
        return Vpst_ruq;
    }
    else{
        return Vpst_ruq;
    }
}

void flexmat6::update_as_pst_uqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vt*iNp*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vr*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pst_uqr(){
    if(Npst_uqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vt*iNs*iNp;
        locations.col(1) = vu + vq*iNu + vr*iNq*iNu;
        Vpst_uqr = sp_mat(locations.t(), vValues, iNp*iNs*iNt,iNu*iNq*iNr);
        Npst_uqr = 1;
        return Vpst_uqr;
    }
    else{
        return Vpst_uqr;
    }
}

void flexmat6::update_as_pst_urq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vt*iNp*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vq*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pst_urq(){
    if(Npst_urq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vt*iNs*iNp;
        locations.col(1) = vu + vr*iNu + vq*iNr*iNu;
        Vpst_urq = sp_mat(locations.t(), vValues, iNp*iNs*iNt,iNu*iNr*iNq);
        Npst_urq = 1;
        return Vpst_urq;
    }
    else{
        return Vpst_urq;
    }
}

void flexmat6::update_as_psu_qrt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vu*iNp*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vt*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psu_qrt(){
    if(Npsu_qrt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vu*iNs*iNp;
        locations.col(1) = vq + vr*iNq + vt*iNr*iNq;
        Vpsu_qrt = sp_mat(locations.t(), vValues, iNp*iNs*iNu,iNq*iNr*iNt);
        Npsu_qrt = 1;
        return Vpsu_qrt;
    }
    else{
        return Vpsu_qrt;
    }
}

void flexmat6::update_as_psu_qtr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vu*iNp*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vr*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psu_qtr(){
    if(Npsu_qtr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vu*iNs*iNp;
        locations.col(1) = vq + vt*iNq + vr*iNt*iNq;
        Vpsu_qtr = sp_mat(locations.t(), vValues, iNp*iNs*iNu,iNq*iNt*iNr);
        Npsu_qtr = 1;
        return Vpsu_qtr;
    }
    else{
        return Vpsu_qtr;
    }
}

void flexmat6::update_as_psu_rqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vu*iNp*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vt*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psu_rqt(){
    if(Npsu_rqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vu*iNs*iNp;
        locations.col(1) = vr + vq*iNr + vt*iNq*iNr;
        Vpsu_rqt = sp_mat(locations.t(), vValues, iNp*iNs*iNu,iNr*iNq*iNt);
        Npsu_rqt = 1;
        return Vpsu_rqt;
    }
    else{
        return Vpsu_rqt;
    }
}

void flexmat6::update_as_psu_rtq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vu*iNp*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vq*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psu_rtq(){
    if(Npsu_rtq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vu*iNs*iNp;
        locations.col(1) = vr + vt*iNr + vq*iNt*iNr;
        Vpsu_rtq = sp_mat(locations.t(), vValues, iNp*iNs*iNu,iNr*iNt*iNq);
        Npsu_rtq = 1;
        return Vpsu_rtq;
    }
    else{
        return Vpsu_rtq;
    }
}

void flexmat6::update_as_psu_tqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vu*iNp*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vr*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psu_tqr(){
    if(Npsu_tqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vu*iNs*iNp;
        locations.col(1) = vt + vq*iNt + vr*iNq*iNt;
        Vpsu_tqr = sp_mat(locations.t(), vValues, iNp*iNs*iNu,iNt*iNq*iNr);
        Npsu_tqr = 1;
        return Vpsu_tqr;
    }
    else{
        return Vpsu_tqr;
    }
}

void flexmat6::update_as_psu_trq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNp - vu*iNp*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vq*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::psu_trq(){
    if(Npsu_trq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vs*iNp + vu*iNs*iNp;
        locations.col(1) = vt + vr*iNt + vq*iNr*iNt;
        Vpsu_trq = sp_mat(locations.t(), vValues, iNp*iNs*iNu,iNt*iNr*iNq);
        Npsu_trq = 1;
        return Vpsu_trq;
    }
    else{
        return Vpsu_trq;
    }
}

void flexmat6::update_as_ptq_rsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vq*iNp*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vu*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptq_rsu(){
    if(Nptq_rsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vq*iNt*iNp;
        locations.col(1) = vr + vs*iNr + vu*iNs*iNr;
        Vptq_rsu = sp_mat(locations.t(), vValues, iNp*iNt*iNq,iNr*iNs*iNu);
        Nptq_rsu = 1;
        return Vptq_rsu;
    }
    else{
        return Vptq_rsu;
    }
}

void flexmat6::update_as_ptq_rus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vq*iNp*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vs*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptq_rus(){
    if(Nptq_rus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vq*iNt*iNp;
        locations.col(1) = vr + vu*iNr + vs*iNu*iNr;
        Vptq_rus = sp_mat(locations.t(), vValues, iNp*iNt*iNq,iNr*iNu*iNs);
        Nptq_rus = 1;
        return Vptq_rus;
    }
    else{
        return Vptq_rus;
    }
}

void flexmat6::update_as_ptq_sru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vq*iNp*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vu*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptq_sru(){
    if(Nptq_sru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vq*iNt*iNp;
        locations.col(1) = vs + vr*iNs + vu*iNr*iNs;
        Vptq_sru = sp_mat(locations.t(), vValues, iNp*iNt*iNq,iNs*iNr*iNu);
        Nptq_sru = 1;
        return Vptq_sru;
    }
    else{
        return Vptq_sru;
    }
}

void flexmat6::update_as_ptq_sur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vq*iNp*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vr*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptq_sur(){
    if(Nptq_sur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vq*iNt*iNp;
        locations.col(1) = vs + vu*iNs + vr*iNu*iNs;
        Vptq_sur = sp_mat(locations.t(), vValues, iNp*iNt*iNq,iNs*iNu*iNr);
        Nptq_sur = 1;
        return Vptq_sur;
    }
    else{
        return Vptq_sur;
    }
}

void flexmat6::update_as_ptq_urs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vq*iNp*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vs*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptq_urs(){
    if(Nptq_urs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vq*iNt*iNp;
        locations.col(1) = vu + vr*iNu + vs*iNr*iNu;
        Vptq_urs = sp_mat(locations.t(), vValues, iNp*iNt*iNq,iNu*iNr*iNs);
        Nptq_urs = 1;
        return Vptq_urs;
    }
    else{
        return Vptq_urs;
    }
}

void flexmat6::update_as_ptq_usr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vq*iNp*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vr*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptq_usr(){
    if(Nptq_usr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vq*iNt*iNp;
        locations.col(1) = vu + vs*iNu + vr*iNs*iNu;
        Vptq_usr = sp_mat(locations.t(), vValues, iNp*iNt*iNq,iNu*iNs*iNr);
        Nptq_usr = 1;
        return Vptq_usr;
    }
    else{
        return Vptq_usr;
    }
}

void flexmat6::update_as_ptr_qsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vr*iNp*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vu*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptr_qsu(){
    if(Nptr_qsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vr*iNt*iNp;
        locations.col(1) = vq + vs*iNq + vu*iNs*iNq;
        Vptr_qsu = sp_mat(locations.t(), vValues, iNp*iNt*iNr,iNq*iNs*iNu);
        Nptr_qsu = 1;
        return Vptr_qsu;
    }
    else{
        return Vptr_qsu;
    }
}

void flexmat6::update_as_ptr_qus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vr*iNp*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vs*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptr_qus(){
    if(Nptr_qus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vr*iNt*iNp;
        locations.col(1) = vq + vu*iNq + vs*iNu*iNq;
        Vptr_qus = sp_mat(locations.t(), vValues, iNp*iNt*iNr,iNq*iNu*iNs);
        Nptr_qus = 1;
        return Vptr_qus;
    }
    else{
        return Vptr_qus;
    }
}

void flexmat6::update_as_ptr_squ(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vr*iNp*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vu*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptr_squ(){
    if(Nptr_squ == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vr*iNt*iNp;
        locations.col(1) = vs + vq*iNs + vu*iNq*iNs;
        Vptr_squ = sp_mat(locations.t(), vValues, iNp*iNt*iNr,iNs*iNq*iNu);
        Nptr_squ = 1;
        return Vptr_squ;
    }
    else{
        return Vptr_squ;
    }
}

void flexmat6::update_as_ptr_suq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vr*iNp*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vq*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptr_suq(){
    if(Nptr_suq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vr*iNt*iNp;
        locations.col(1) = vs + vu*iNs + vq*iNu*iNs;
        Vptr_suq = sp_mat(locations.t(), vValues, iNp*iNt*iNr,iNs*iNu*iNq);
        Nptr_suq = 1;
        return Vptr_suq;
    }
    else{
        return Vptr_suq;
    }
}

void flexmat6::update_as_ptr_uqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vr*iNp*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vs*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptr_uqs(){
    if(Nptr_uqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vr*iNt*iNp;
        locations.col(1) = vu + vq*iNu + vs*iNq*iNu;
        Vptr_uqs = sp_mat(locations.t(), vValues, iNp*iNt*iNr,iNu*iNq*iNs);
        Nptr_uqs = 1;
        return Vptr_uqs;
    }
    else{
        return Vptr_uqs;
    }
}

void flexmat6::update_as_ptr_usq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vr*iNp*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vq*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptr_usq(){
    if(Nptr_usq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vr*iNt*iNp;
        locations.col(1) = vu + vs*iNu + vq*iNs*iNu;
        Vptr_usq = sp_mat(locations.t(), vValues, iNp*iNt*iNr,iNu*iNs*iNq);
        Nptr_usq = 1;
        return Vptr_usq;
    }
    else{
        return Vptr_usq;
    }
}

void flexmat6::update_as_pts_qru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vs*iNp*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vu*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pts_qru(){
    if(Npts_qru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vs*iNt*iNp;
        locations.col(1) = vq + vr*iNq + vu*iNr*iNq;
        Vpts_qru = sp_mat(locations.t(), vValues, iNp*iNt*iNs,iNq*iNr*iNu);
        Npts_qru = 1;
        return Vpts_qru;
    }
    else{
        return Vpts_qru;
    }
}

void flexmat6::update_as_pts_qur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vs*iNp*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vr*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pts_qur(){
    if(Npts_qur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vs*iNt*iNp;
        locations.col(1) = vq + vu*iNq + vr*iNu*iNq;
        Vpts_qur = sp_mat(locations.t(), vValues, iNp*iNt*iNs,iNq*iNu*iNr);
        Npts_qur = 1;
        return Vpts_qur;
    }
    else{
        return Vpts_qur;
    }
}

void flexmat6::update_as_pts_rqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vs*iNp*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vu*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pts_rqu(){
    if(Npts_rqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vs*iNt*iNp;
        locations.col(1) = vr + vq*iNr + vu*iNq*iNr;
        Vpts_rqu = sp_mat(locations.t(), vValues, iNp*iNt*iNs,iNr*iNq*iNu);
        Npts_rqu = 1;
        return Vpts_rqu;
    }
    else{
        return Vpts_rqu;
    }
}

void flexmat6::update_as_pts_ruq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vs*iNp*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vq*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pts_ruq(){
    if(Npts_ruq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vs*iNt*iNp;
        locations.col(1) = vr + vu*iNr + vq*iNu*iNr;
        Vpts_ruq = sp_mat(locations.t(), vValues, iNp*iNt*iNs,iNr*iNu*iNq);
        Npts_ruq = 1;
        return Vpts_ruq;
    }
    else{
        return Vpts_ruq;
    }
}

void flexmat6::update_as_pts_uqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vs*iNp*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vr*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pts_uqr(){
    if(Npts_uqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vs*iNt*iNp;
        locations.col(1) = vu + vq*iNu + vr*iNq*iNu;
        Vpts_uqr = sp_mat(locations.t(), vValues, iNp*iNt*iNs,iNu*iNq*iNr);
        Npts_uqr = 1;
        return Vpts_uqr;
    }
    else{
        return Vpts_uqr;
    }
}

void flexmat6::update_as_pts_urq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vs*iNp*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vq*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pts_urq(){
    if(Npts_urq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vs*iNt*iNp;
        locations.col(1) = vu + vr*iNu + vq*iNr*iNu;
        Vpts_urq = sp_mat(locations.t(), vValues, iNp*iNt*iNs,iNu*iNr*iNq);
        Npts_urq = 1;
        return Vpts_urq;
    }
    else{
        return Vpts_urq;
    }
}

void flexmat6::update_as_ptu_qrs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vu*iNp*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vs*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptu_qrs(){
    if(Nptu_qrs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vu*iNt*iNp;
        locations.col(1) = vq + vr*iNq + vs*iNr*iNq;
        Vptu_qrs = sp_mat(locations.t(), vValues, iNp*iNt*iNu,iNq*iNr*iNs);
        Nptu_qrs = 1;
        return Vptu_qrs;
    }
    else{
        return Vptu_qrs;
    }
}

void flexmat6::update_as_ptu_qsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vu*iNp*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vr*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptu_qsr(){
    if(Nptu_qsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vu*iNt*iNp;
        locations.col(1) = vq + vs*iNq + vr*iNs*iNq;
        Vptu_qsr = sp_mat(locations.t(), vValues, iNp*iNt*iNu,iNq*iNs*iNr);
        Nptu_qsr = 1;
        return Vptu_qsr;
    }
    else{
        return Vptu_qsr;
    }
}

void flexmat6::update_as_ptu_rqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vu*iNp*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vs*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptu_rqs(){
    if(Nptu_rqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vu*iNt*iNp;
        locations.col(1) = vr + vq*iNr + vs*iNq*iNr;
        Vptu_rqs = sp_mat(locations.t(), vValues, iNp*iNt*iNu,iNr*iNq*iNs);
        Nptu_rqs = 1;
        return Vptu_rqs;
    }
    else{
        return Vptu_rqs;
    }
}

void flexmat6::update_as_ptu_rsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vu*iNp*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vq*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptu_rsq(){
    if(Nptu_rsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vu*iNt*iNp;
        locations.col(1) = vr + vs*iNr + vq*iNs*iNr;
        Vptu_rsq = sp_mat(locations.t(), vValues, iNp*iNt*iNu,iNr*iNs*iNq);
        Nptu_rsq = 1;
        return Vptu_rsq;
    }
    else{
        return Vptu_rsq;
    }
}

void flexmat6::update_as_ptu_sqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vu*iNp*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vr*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptu_sqr(){
    if(Nptu_sqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vu*iNt*iNp;
        locations.col(1) = vs + vq*iNs + vr*iNq*iNs;
        Vptu_sqr = sp_mat(locations.t(), vValues, iNp*iNt*iNu,iNs*iNq*iNr);
        Nptu_sqr = 1;
        return Vptu_sqr;
    }
    else{
        return Vptu_sqr;
    }
}

void flexmat6::update_as_ptu_srq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNp - vu*iNp*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vq*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ptu_srq(){
    if(Nptu_srq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vt*iNp + vu*iNt*iNp;
        locations.col(1) = vs + vr*iNs + vq*iNr*iNs;
        Vptu_srq = sp_mat(locations.t(), vValues, iNp*iNt*iNu,iNs*iNr*iNq);
        Nptu_srq = 1;
        return Vptu_srq;
    }
    else{
        return Vptu_srq;
    }
}

void flexmat6::update_as_puq_rst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vq*iNp*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vt*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::puq_rst(){
    if(Npuq_rst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vq*iNu*iNp;
        locations.col(1) = vr + vs*iNr + vt*iNs*iNr;
        Vpuq_rst = sp_mat(locations.t(), vValues, iNp*iNu*iNq,iNr*iNs*iNt);
        Npuq_rst = 1;
        return Vpuq_rst;
    }
    else{
        return Vpuq_rst;
    }
}

void flexmat6::update_as_puq_rts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vq*iNp*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vs*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::puq_rts(){
    if(Npuq_rts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vq*iNu*iNp;
        locations.col(1) = vr + vt*iNr + vs*iNt*iNr;
        Vpuq_rts = sp_mat(locations.t(), vValues, iNp*iNu*iNq,iNr*iNt*iNs);
        Npuq_rts = 1;
        return Vpuq_rts;
    }
    else{
        return Vpuq_rts;
    }
}

void flexmat6::update_as_puq_srt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vq*iNp*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vt*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::puq_srt(){
    if(Npuq_srt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vq*iNu*iNp;
        locations.col(1) = vs + vr*iNs + vt*iNr*iNs;
        Vpuq_srt = sp_mat(locations.t(), vValues, iNp*iNu*iNq,iNs*iNr*iNt);
        Npuq_srt = 1;
        return Vpuq_srt;
    }
    else{
        return Vpuq_srt;
    }
}

void flexmat6::update_as_puq_str(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vq*iNp*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vr*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::puq_str(){
    if(Npuq_str == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vq*iNu*iNp;
        locations.col(1) = vs + vt*iNs + vr*iNt*iNs;
        Vpuq_str = sp_mat(locations.t(), vValues, iNp*iNu*iNq,iNs*iNt*iNr);
        Npuq_str = 1;
        return Vpuq_str;
    }
    else{
        return Vpuq_str;
    }
}

void flexmat6::update_as_puq_trs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vq*iNp*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vs*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::puq_trs(){
    if(Npuq_trs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vq*iNu*iNp;
        locations.col(1) = vt + vr*iNt + vs*iNr*iNt;
        Vpuq_trs = sp_mat(locations.t(), vValues, iNp*iNu*iNq,iNt*iNr*iNs);
        Npuq_trs = 1;
        return Vpuq_trs;
    }
    else{
        return Vpuq_trs;
    }
}

void flexmat6::update_as_puq_tsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vq*iNp*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vr*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::puq_tsr(){
    if(Npuq_tsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vq*iNu*iNp;
        locations.col(1) = vt + vs*iNt + vr*iNs*iNt;
        Vpuq_tsr = sp_mat(locations.t(), vValues, iNp*iNu*iNq,iNt*iNs*iNr);
        Npuq_tsr = 1;
        return Vpuq_tsr;
    }
    else{
        return Vpuq_tsr;
    }
}

void flexmat6::update_as_pur_qst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vr*iNp*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vt*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pur_qst(){
    if(Npur_qst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vr*iNu*iNp;
        locations.col(1) = vq + vs*iNq + vt*iNs*iNq;
        Vpur_qst = sp_mat(locations.t(), vValues, iNp*iNu*iNr,iNq*iNs*iNt);
        Npur_qst = 1;
        return Vpur_qst;
    }
    else{
        return Vpur_qst;
    }
}

void flexmat6::update_as_pur_qts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vr*iNp*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vs*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pur_qts(){
    if(Npur_qts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vr*iNu*iNp;
        locations.col(1) = vq + vt*iNq + vs*iNt*iNq;
        Vpur_qts = sp_mat(locations.t(), vValues, iNp*iNu*iNr,iNq*iNt*iNs);
        Npur_qts = 1;
        return Vpur_qts;
    }
    else{
        return Vpur_qts;
    }
}

void flexmat6::update_as_pur_sqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vr*iNp*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vt*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pur_sqt(){
    if(Npur_sqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vr*iNu*iNp;
        locations.col(1) = vs + vq*iNs + vt*iNq*iNs;
        Vpur_sqt = sp_mat(locations.t(), vValues, iNp*iNu*iNr,iNs*iNq*iNt);
        Npur_sqt = 1;
        return Vpur_sqt;
    }
    else{
        return Vpur_sqt;
    }
}

void flexmat6::update_as_pur_stq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vr*iNp*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vq*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pur_stq(){
    if(Npur_stq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vr*iNu*iNp;
        locations.col(1) = vs + vt*iNs + vq*iNt*iNs;
        Vpur_stq = sp_mat(locations.t(), vValues, iNp*iNu*iNr,iNs*iNt*iNq);
        Npur_stq = 1;
        return Vpur_stq;
    }
    else{
        return Vpur_stq;
    }
}

void flexmat6::update_as_pur_tqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vr*iNp*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vs*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pur_tqs(){
    if(Npur_tqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vr*iNu*iNp;
        locations.col(1) = vt + vq*iNt + vs*iNq*iNt;
        Vpur_tqs = sp_mat(locations.t(), vValues, iNp*iNu*iNr,iNt*iNq*iNs);
        Npur_tqs = 1;
        return Vpur_tqs;
    }
    else{
        return Vpur_tqs;
    }
}

void flexmat6::update_as_pur_tsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vr*iNp*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vq*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pur_tsq(){
    if(Npur_tsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vr*iNu*iNp;
        locations.col(1) = vt + vs*iNt + vq*iNs*iNt;
        Vpur_tsq = sp_mat(locations.t(), vValues, iNp*iNu*iNr,iNt*iNs*iNq);
        Npur_tsq = 1;
        return Vpur_tsq;
    }
    else{
        return Vpur_tsq;
    }
}

void flexmat6::update_as_pus_qrt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vs*iNp*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vt*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pus_qrt(){
    if(Npus_qrt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vs*iNu*iNp;
        locations.col(1) = vq + vr*iNq + vt*iNr*iNq;
        Vpus_qrt = sp_mat(locations.t(), vValues, iNp*iNu*iNs,iNq*iNr*iNt);
        Npus_qrt = 1;
        return Vpus_qrt;
    }
    else{
        return Vpus_qrt;
    }
}

void flexmat6::update_as_pus_qtr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vs*iNp*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vr*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pus_qtr(){
    if(Npus_qtr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vs*iNu*iNp;
        locations.col(1) = vq + vt*iNq + vr*iNt*iNq;
        Vpus_qtr = sp_mat(locations.t(), vValues, iNp*iNu*iNs,iNq*iNt*iNr);
        Npus_qtr = 1;
        return Vpus_qtr;
    }
    else{
        return Vpus_qtr;
    }
}

void flexmat6::update_as_pus_rqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vs*iNp*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vt*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pus_rqt(){
    if(Npus_rqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vs*iNu*iNp;
        locations.col(1) = vr + vq*iNr + vt*iNq*iNr;
        Vpus_rqt = sp_mat(locations.t(), vValues, iNp*iNu*iNs,iNr*iNq*iNt);
        Npus_rqt = 1;
        return Vpus_rqt;
    }
    else{
        return Vpus_rqt;
    }
}

void flexmat6::update_as_pus_rtq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vs*iNp*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vq*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pus_rtq(){
    if(Npus_rtq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vs*iNu*iNp;
        locations.col(1) = vr + vt*iNr + vq*iNt*iNr;
        Vpus_rtq = sp_mat(locations.t(), vValues, iNp*iNu*iNs,iNr*iNt*iNq);
        Npus_rtq = 1;
        return Vpus_rtq;
    }
    else{
        return Vpus_rtq;
    }
}

void flexmat6::update_as_pus_tqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vs*iNp*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vr*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pus_tqr(){
    if(Npus_tqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vs*iNu*iNp;
        locations.col(1) = vt + vq*iNt + vr*iNq*iNt;
        Vpus_tqr = sp_mat(locations.t(), vValues, iNp*iNu*iNs,iNt*iNq*iNr);
        Npus_tqr = 1;
        return Vpus_tqr;
    }
    else{
        return Vpus_tqr;
    }
}

void flexmat6::update_as_pus_trq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vs*iNp*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vq*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::pus_trq(){
    if(Npus_trq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vs*iNu*iNp;
        locations.col(1) = vt + vr*iNt + vq*iNr*iNt;
        Vpus_trq = sp_mat(locations.t(), vValues, iNp*iNu*iNs,iNt*iNr*iNq);
        Npus_trq = 1;
        return Vpus_trq;
    }
    else{
        return Vpus_trq;
    }
}

void flexmat6::update_as_put_qrs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vt*iNp*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vs*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::put_qrs(){
    if(Nput_qrs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vt*iNu*iNp;
        locations.col(1) = vq + vr*iNq + vs*iNr*iNq;
        Vput_qrs = sp_mat(locations.t(), vValues, iNp*iNu*iNt,iNq*iNr*iNs);
        Nput_qrs = 1;
        return Vput_qrs;
    }
    else{
        return Vput_qrs;
    }
}

void flexmat6::update_as_put_qsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vt*iNp*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vr*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::put_qsr(){
    if(Nput_qsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vt*iNu*iNp;
        locations.col(1) = vq + vs*iNq + vr*iNs*iNq;
        Vput_qsr = sp_mat(locations.t(), vValues, iNp*iNu*iNt,iNq*iNs*iNr);
        Nput_qsr = 1;
        return Vput_qsr;
    }
    else{
        return Vput_qsr;
    }
}

void flexmat6::update_as_put_rqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vt*iNp*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vs*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::put_rqs(){
    if(Nput_rqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vt*iNu*iNp;
        locations.col(1) = vr + vq*iNr + vs*iNq*iNr;
        Vput_rqs = sp_mat(locations.t(), vValues, iNp*iNu*iNt,iNr*iNq*iNs);
        Nput_rqs = 1;
        return Vput_rqs;
    }
    else{
        return Vput_rqs;
    }
}

void flexmat6::update_as_put_rsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vt*iNp*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vq*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::put_rsq(){
    if(Nput_rsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vt*iNu*iNp;
        locations.col(1) = vr + vs*iNr + vq*iNs*iNr;
        Vput_rsq = sp_mat(locations.t(), vValues, iNp*iNu*iNt,iNr*iNs*iNq);
        Nput_rsq = 1;
        return Vput_rsq;
    }
    else{
        return Vput_rsq;
    }
}

void flexmat6::update_as_put_sqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vt*iNp*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vr*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::put_sqr(){
    if(Nput_sqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vt*iNu*iNp;
        locations.col(1) = vs + vq*iNs + vr*iNq*iNs;
        Vput_sqr = sp_mat(locations.t(), vValues, iNp*iNu*iNt,iNs*iNq*iNr);
        Nput_sqr = 1;
        return Vput_sqr;
    }
    else{
        return Vput_sqr;
    }
}

void flexmat6::update_as_put_srq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNp - vt*iNp*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vq*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::put_srq(){
    if(Nput_srq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vp + vu*iNp + vt*iNu*iNp;
        locations.col(1) = vs + vr*iNs + vq*iNr*iNs;
        Vput_srq = sp_mat(locations.t(), vValues, iNp*iNu*iNt,iNs*iNr*iNq);
        Nput_srq = 1;
        return Vput_srq;
    }
    else{
        return Vput_srq;
    }
}

void flexmat6::update_as_qpr_stu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vr*iNq*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vu*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpr_stu(){
    if(Nqpr_stu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vr*iNp*iNq;
        locations.col(1) = vs + vt*iNs + vu*iNt*iNs;
        Vqpr_stu = sp_mat(locations.t(), vValues, iNq*iNp*iNr,iNs*iNt*iNu);
        Nqpr_stu = 1;
        return Vqpr_stu;
    }
    else{
        return Vqpr_stu;
    }
}

void flexmat6::update_as_qpr_sut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vr*iNq*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vt*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpr_sut(){
    if(Nqpr_sut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vr*iNp*iNq;
        locations.col(1) = vs + vu*iNs + vt*iNu*iNs;
        Vqpr_sut = sp_mat(locations.t(), vValues, iNq*iNp*iNr,iNs*iNu*iNt);
        Nqpr_sut = 1;
        return Vqpr_sut;
    }
    else{
        return Vqpr_sut;
    }
}

void flexmat6::update_as_qpr_tsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vr*iNq*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vu*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpr_tsu(){
    if(Nqpr_tsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vr*iNp*iNq;
        locations.col(1) = vt + vs*iNt + vu*iNs*iNt;
        Vqpr_tsu = sp_mat(locations.t(), vValues, iNq*iNp*iNr,iNt*iNs*iNu);
        Nqpr_tsu = 1;
        return Vqpr_tsu;
    }
    else{
        return Vqpr_tsu;
    }
}

void flexmat6::update_as_qpr_tus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vr*iNq*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vs*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpr_tus(){
    if(Nqpr_tus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vr*iNp*iNq;
        locations.col(1) = vt + vu*iNt + vs*iNu*iNt;
        Vqpr_tus = sp_mat(locations.t(), vValues, iNq*iNp*iNr,iNt*iNu*iNs);
        Nqpr_tus = 1;
        return Vqpr_tus;
    }
    else{
        return Vqpr_tus;
    }
}

void flexmat6::update_as_qpr_ust(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vr*iNq*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vt*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpr_ust(){
    if(Nqpr_ust == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vr*iNp*iNq;
        locations.col(1) = vu + vs*iNu + vt*iNs*iNu;
        Vqpr_ust = sp_mat(locations.t(), vValues, iNq*iNp*iNr,iNu*iNs*iNt);
        Nqpr_ust = 1;
        return Vqpr_ust;
    }
    else{
        return Vqpr_ust;
    }
}

void flexmat6::update_as_qpr_uts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vr*iNq*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vs*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpr_uts(){
    if(Nqpr_uts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vr*iNp*iNq;
        locations.col(1) = vu + vt*iNu + vs*iNt*iNu;
        Vqpr_uts = sp_mat(locations.t(), vValues, iNq*iNp*iNr,iNu*iNt*iNs);
        Nqpr_uts = 1;
        return Vqpr_uts;
    }
    else{
        return Vqpr_uts;
    }
}

void flexmat6::update_as_qps_rtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vs*iNq*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vu*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qps_rtu(){
    if(Nqps_rtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vs*iNp*iNq;
        locations.col(1) = vr + vt*iNr + vu*iNt*iNr;
        Vqps_rtu = sp_mat(locations.t(), vValues, iNq*iNp*iNs,iNr*iNt*iNu);
        Nqps_rtu = 1;
        return Vqps_rtu;
    }
    else{
        return Vqps_rtu;
    }
}

void flexmat6::update_as_qps_rut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vs*iNq*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vt*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qps_rut(){
    if(Nqps_rut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vs*iNp*iNq;
        locations.col(1) = vr + vu*iNr + vt*iNu*iNr;
        Vqps_rut = sp_mat(locations.t(), vValues, iNq*iNp*iNs,iNr*iNu*iNt);
        Nqps_rut = 1;
        return Vqps_rut;
    }
    else{
        return Vqps_rut;
    }
}

void flexmat6::update_as_qps_tru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vs*iNq*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vu*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qps_tru(){
    if(Nqps_tru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vs*iNp*iNq;
        locations.col(1) = vt + vr*iNt + vu*iNr*iNt;
        Vqps_tru = sp_mat(locations.t(), vValues, iNq*iNp*iNs,iNt*iNr*iNu);
        Nqps_tru = 1;
        return Vqps_tru;
    }
    else{
        return Vqps_tru;
    }
}

void flexmat6::update_as_qps_tur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vs*iNq*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vr*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qps_tur(){
    if(Nqps_tur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vs*iNp*iNq;
        locations.col(1) = vt + vu*iNt + vr*iNu*iNt;
        Vqps_tur = sp_mat(locations.t(), vValues, iNq*iNp*iNs,iNt*iNu*iNr);
        Nqps_tur = 1;
        return Vqps_tur;
    }
    else{
        return Vqps_tur;
    }
}

void flexmat6::update_as_qps_urt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vs*iNq*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vt*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qps_urt(){
    if(Nqps_urt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vs*iNp*iNq;
        locations.col(1) = vu + vr*iNu + vt*iNr*iNu;
        Vqps_urt = sp_mat(locations.t(), vValues, iNq*iNp*iNs,iNu*iNr*iNt);
        Nqps_urt = 1;
        return Vqps_urt;
    }
    else{
        return Vqps_urt;
    }
}

void flexmat6::update_as_qps_utr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vs*iNq*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vr*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qps_utr(){
    if(Nqps_utr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vs*iNp*iNq;
        locations.col(1) = vu + vt*iNu + vr*iNt*iNu;
        Vqps_utr = sp_mat(locations.t(), vValues, iNq*iNp*iNs,iNu*iNt*iNr);
        Nqps_utr = 1;
        return Vqps_utr;
    }
    else{
        return Vqps_utr;
    }
}

void flexmat6::update_as_qpt_rsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vt*iNq*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vu*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpt_rsu(){
    if(Nqpt_rsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vt*iNp*iNq;
        locations.col(1) = vr + vs*iNr + vu*iNs*iNr;
        Vqpt_rsu = sp_mat(locations.t(), vValues, iNq*iNp*iNt,iNr*iNs*iNu);
        Nqpt_rsu = 1;
        return Vqpt_rsu;
    }
    else{
        return Vqpt_rsu;
    }
}

void flexmat6::update_as_qpt_rus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vt*iNq*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vs*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpt_rus(){
    if(Nqpt_rus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vt*iNp*iNq;
        locations.col(1) = vr + vu*iNr + vs*iNu*iNr;
        Vqpt_rus = sp_mat(locations.t(), vValues, iNq*iNp*iNt,iNr*iNu*iNs);
        Nqpt_rus = 1;
        return Vqpt_rus;
    }
    else{
        return Vqpt_rus;
    }
}

void flexmat6::update_as_qpt_sru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vt*iNq*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vu*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpt_sru(){
    if(Nqpt_sru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vt*iNp*iNq;
        locations.col(1) = vs + vr*iNs + vu*iNr*iNs;
        Vqpt_sru = sp_mat(locations.t(), vValues, iNq*iNp*iNt,iNs*iNr*iNu);
        Nqpt_sru = 1;
        return Vqpt_sru;
    }
    else{
        return Vqpt_sru;
    }
}

void flexmat6::update_as_qpt_sur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vt*iNq*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vr*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpt_sur(){
    if(Nqpt_sur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vt*iNp*iNq;
        locations.col(1) = vs + vu*iNs + vr*iNu*iNs;
        Vqpt_sur = sp_mat(locations.t(), vValues, iNq*iNp*iNt,iNs*iNu*iNr);
        Nqpt_sur = 1;
        return Vqpt_sur;
    }
    else{
        return Vqpt_sur;
    }
}

void flexmat6::update_as_qpt_urs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vt*iNq*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vs*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpt_urs(){
    if(Nqpt_urs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vt*iNp*iNq;
        locations.col(1) = vu + vr*iNu + vs*iNr*iNu;
        Vqpt_urs = sp_mat(locations.t(), vValues, iNq*iNp*iNt,iNu*iNr*iNs);
        Nqpt_urs = 1;
        return Vqpt_urs;
    }
    else{
        return Vqpt_urs;
    }
}

void flexmat6::update_as_qpt_usr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vt*iNq*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vr*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpt_usr(){
    if(Nqpt_usr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vt*iNp*iNq;
        locations.col(1) = vu + vs*iNu + vr*iNs*iNu;
        Vqpt_usr = sp_mat(locations.t(), vValues, iNq*iNp*iNt,iNu*iNs*iNr);
        Nqpt_usr = 1;
        return Vqpt_usr;
    }
    else{
        return Vqpt_usr;
    }
}

void flexmat6::update_as_qpu_rst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vu*iNq*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vt*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpu_rst(){
    if(Nqpu_rst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vu*iNp*iNq;
        locations.col(1) = vr + vs*iNr + vt*iNs*iNr;
        Vqpu_rst = sp_mat(locations.t(), vValues, iNq*iNp*iNu,iNr*iNs*iNt);
        Nqpu_rst = 1;
        return Vqpu_rst;
    }
    else{
        return Vqpu_rst;
    }
}

void flexmat6::update_as_qpu_rts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vu*iNq*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vs*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpu_rts(){
    if(Nqpu_rts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vu*iNp*iNq;
        locations.col(1) = vr + vt*iNr + vs*iNt*iNr;
        Vqpu_rts = sp_mat(locations.t(), vValues, iNq*iNp*iNu,iNr*iNt*iNs);
        Nqpu_rts = 1;
        return Vqpu_rts;
    }
    else{
        return Vqpu_rts;
    }
}

void flexmat6::update_as_qpu_srt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vu*iNq*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vt*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpu_srt(){
    if(Nqpu_srt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vu*iNp*iNq;
        locations.col(1) = vs + vr*iNs + vt*iNr*iNs;
        Vqpu_srt = sp_mat(locations.t(), vValues, iNq*iNp*iNu,iNs*iNr*iNt);
        Nqpu_srt = 1;
        return Vqpu_srt;
    }
    else{
        return Vqpu_srt;
    }
}

void flexmat6::update_as_qpu_str(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vu*iNq*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vr*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpu_str(){
    if(Nqpu_str == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vu*iNp*iNq;
        locations.col(1) = vs + vt*iNs + vr*iNt*iNs;
        Vqpu_str = sp_mat(locations.t(), vValues, iNq*iNp*iNu,iNs*iNt*iNr);
        Nqpu_str = 1;
        return Vqpu_str;
    }
    else{
        return Vqpu_str;
    }
}

void flexmat6::update_as_qpu_trs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vu*iNq*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vs*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpu_trs(){
    if(Nqpu_trs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vu*iNp*iNq;
        locations.col(1) = vt + vr*iNt + vs*iNr*iNt;
        Vqpu_trs = sp_mat(locations.t(), vValues, iNq*iNp*iNu,iNt*iNr*iNs);
        Nqpu_trs = 1;
        return Vqpu_trs;
    }
    else{
        return Vqpu_trs;
    }
}

void flexmat6::update_as_qpu_tsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNq - vu*iNq*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vr*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qpu_tsr(){
    if(Nqpu_tsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vp*iNq + vu*iNp*iNq;
        locations.col(1) = vt + vs*iNt + vr*iNs*iNt;
        Vqpu_tsr = sp_mat(locations.t(), vValues, iNq*iNp*iNu,iNt*iNs*iNr);
        Nqpu_tsr = 1;
        return Vqpu_tsr;
    }
    else{
        return Vqpu_tsr;
    }
}

void flexmat6::update_as_qrp_stu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vp*iNq*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vu*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrp_stu(){
    if(Nqrp_stu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vp*iNr*iNq;
        locations.col(1) = vs + vt*iNs + vu*iNt*iNs;
        Vqrp_stu = sp_mat(locations.t(), vValues, iNq*iNr*iNp,iNs*iNt*iNu);
        Nqrp_stu = 1;
        return Vqrp_stu;
    }
    else{
        return Vqrp_stu;
    }
}

void flexmat6::update_as_qrp_sut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vp*iNq*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vt*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrp_sut(){
    if(Nqrp_sut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vp*iNr*iNq;
        locations.col(1) = vs + vu*iNs + vt*iNu*iNs;
        Vqrp_sut = sp_mat(locations.t(), vValues, iNq*iNr*iNp,iNs*iNu*iNt);
        Nqrp_sut = 1;
        return Vqrp_sut;
    }
    else{
        return Vqrp_sut;
    }
}

void flexmat6::update_as_qrp_tsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vp*iNq*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vu*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrp_tsu(){
    if(Nqrp_tsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vp*iNr*iNq;
        locations.col(1) = vt + vs*iNt + vu*iNs*iNt;
        Vqrp_tsu = sp_mat(locations.t(), vValues, iNq*iNr*iNp,iNt*iNs*iNu);
        Nqrp_tsu = 1;
        return Vqrp_tsu;
    }
    else{
        return Vqrp_tsu;
    }
}

void flexmat6::update_as_qrp_tus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vp*iNq*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vs*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrp_tus(){
    if(Nqrp_tus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vp*iNr*iNq;
        locations.col(1) = vt + vu*iNt + vs*iNu*iNt;
        Vqrp_tus = sp_mat(locations.t(), vValues, iNq*iNr*iNp,iNt*iNu*iNs);
        Nqrp_tus = 1;
        return Vqrp_tus;
    }
    else{
        return Vqrp_tus;
    }
}

void flexmat6::update_as_qrp_ust(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vp*iNq*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vt*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrp_ust(){
    if(Nqrp_ust == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vp*iNr*iNq;
        locations.col(1) = vu + vs*iNu + vt*iNs*iNu;
        Vqrp_ust = sp_mat(locations.t(), vValues, iNq*iNr*iNp,iNu*iNs*iNt);
        Nqrp_ust = 1;
        return Vqrp_ust;
    }
    else{
        return Vqrp_ust;
    }
}

void flexmat6::update_as_qrp_uts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vp*iNq*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vs*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrp_uts(){
    if(Nqrp_uts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vp*iNr*iNq;
        locations.col(1) = vu + vt*iNu + vs*iNt*iNu;
        Vqrp_uts = sp_mat(locations.t(), vValues, iNq*iNr*iNp,iNu*iNt*iNs);
        Nqrp_uts = 1;
        return Vqrp_uts;
    }
    else{
        return Vqrp_uts;
    }
}

void flexmat6::update_as_qrs_ptu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vs*iNq*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vu*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrs_ptu(){
    if(Nqrs_ptu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vs*iNr*iNq;
        locations.col(1) = vp + vt*iNp + vu*iNt*iNp;
        Vqrs_ptu = sp_mat(locations.t(), vValues, iNq*iNr*iNs,iNp*iNt*iNu);
        Nqrs_ptu = 1;
        return Vqrs_ptu;
    }
    else{
        return Vqrs_ptu;
    }
}

void flexmat6::update_as_qrs_put(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vs*iNq*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vt*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrs_put(){
    if(Nqrs_put == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vs*iNr*iNq;
        locations.col(1) = vp + vu*iNp + vt*iNu*iNp;
        Vqrs_put = sp_mat(locations.t(), vValues, iNq*iNr*iNs,iNp*iNu*iNt);
        Nqrs_put = 1;
        return Vqrs_put;
    }
    else{
        return Vqrs_put;
    }
}

void flexmat6::update_as_qrs_tpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vs*iNq*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vu*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrs_tpu(){
    if(Nqrs_tpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vs*iNr*iNq;
        locations.col(1) = vt + vp*iNt + vu*iNp*iNt;
        Vqrs_tpu = sp_mat(locations.t(), vValues, iNq*iNr*iNs,iNt*iNp*iNu);
        Nqrs_tpu = 1;
        return Vqrs_tpu;
    }
    else{
        return Vqrs_tpu;
    }
}

void flexmat6::update_as_qrs_tup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vs*iNq*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vp*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrs_tup(){
    if(Nqrs_tup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vs*iNr*iNq;
        locations.col(1) = vt + vu*iNt + vp*iNu*iNt;
        Vqrs_tup = sp_mat(locations.t(), vValues, iNq*iNr*iNs,iNt*iNu*iNp);
        Nqrs_tup = 1;
        return Vqrs_tup;
    }
    else{
        return Vqrs_tup;
    }
}

void flexmat6::update_as_qrs_upt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vs*iNq*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vt*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrs_upt(){
    if(Nqrs_upt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vs*iNr*iNq;
        locations.col(1) = vu + vp*iNu + vt*iNp*iNu;
        Vqrs_upt = sp_mat(locations.t(), vValues, iNq*iNr*iNs,iNu*iNp*iNt);
        Nqrs_upt = 1;
        return Vqrs_upt;
    }
    else{
        return Vqrs_upt;
    }
}

void flexmat6::update_as_qrs_utp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vs*iNq*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vp*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrs_utp(){
    if(Nqrs_utp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vs*iNr*iNq;
        locations.col(1) = vu + vt*iNu + vp*iNt*iNu;
        Vqrs_utp = sp_mat(locations.t(), vValues, iNq*iNr*iNs,iNu*iNt*iNp);
        Nqrs_utp = 1;
        return Vqrs_utp;
    }
    else{
        return Vqrs_utp;
    }
}

void flexmat6::update_as_qrt_psu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vt*iNq*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vu*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrt_psu(){
    if(Nqrt_psu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vt*iNr*iNq;
        locations.col(1) = vp + vs*iNp + vu*iNs*iNp;
        Vqrt_psu = sp_mat(locations.t(), vValues, iNq*iNr*iNt,iNp*iNs*iNu);
        Nqrt_psu = 1;
        return Vqrt_psu;
    }
    else{
        return Vqrt_psu;
    }
}

void flexmat6::update_as_qrt_pus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vt*iNq*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vs*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrt_pus(){
    if(Nqrt_pus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vt*iNr*iNq;
        locations.col(1) = vp + vu*iNp + vs*iNu*iNp;
        Vqrt_pus = sp_mat(locations.t(), vValues, iNq*iNr*iNt,iNp*iNu*iNs);
        Nqrt_pus = 1;
        return Vqrt_pus;
    }
    else{
        return Vqrt_pus;
    }
}

void flexmat6::update_as_qrt_spu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vt*iNq*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vu*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrt_spu(){
    if(Nqrt_spu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vt*iNr*iNq;
        locations.col(1) = vs + vp*iNs + vu*iNp*iNs;
        Vqrt_spu = sp_mat(locations.t(), vValues, iNq*iNr*iNt,iNs*iNp*iNu);
        Nqrt_spu = 1;
        return Vqrt_spu;
    }
    else{
        return Vqrt_spu;
    }
}

void flexmat6::update_as_qrt_sup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vt*iNq*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vp*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrt_sup(){
    if(Nqrt_sup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vt*iNr*iNq;
        locations.col(1) = vs + vu*iNs + vp*iNu*iNs;
        Vqrt_sup = sp_mat(locations.t(), vValues, iNq*iNr*iNt,iNs*iNu*iNp);
        Nqrt_sup = 1;
        return Vqrt_sup;
    }
    else{
        return Vqrt_sup;
    }
}

void flexmat6::update_as_qrt_ups(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vt*iNq*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vs*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrt_ups(){
    if(Nqrt_ups == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vt*iNr*iNq;
        locations.col(1) = vu + vp*iNu + vs*iNp*iNu;
        Vqrt_ups = sp_mat(locations.t(), vValues, iNq*iNr*iNt,iNu*iNp*iNs);
        Nqrt_ups = 1;
        return Vqrt_ups;
    }
    else{
        return Vqrt_ups;
    }
}

void flexmat6::update_as_qrt_usp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vt*iNq*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vp*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qrt_usp(){
    if(Nqrt_usp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vt*iNr*iNq;
        locations.col(1) = vu + vs*iNu + vp*iNs*iNu;
        Vqrt_usp = sp_mat(locations.t(), vValues, iNq*iNr*iNt,iNu*iNs*iNp);
        Nqrt_usp = 1;
        return Vqrt_usp;
    }
    else{
        return Vqrt_usp;
    }
}

void flexmat6::update_as_qru_pst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vu*iNq*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vt*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qru_pst(){
    if(Nqru_pst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vu*iNr*iNq;
        locations.col(1) = vp + vs*iNp + vt*iNs*iNp;
        Vqru_pst = sp_mat(locations.t(), vValues, iNq*iNr*iNu,iNp*iNs*iNt);
        Nqru_pst = 1;
        return Vqru_pst;
    }
    else{
        return Vqru_pst;
    }
}

void flexmat6::update_as_qru_pts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vu*iNq*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vs*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qru_pts(){
    if(Nqru_pts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vu*iNr*iNq;
        locations.col(1) = vp + vt*iNp + vs*iNt*iNp;
        Vqru_pts = sp_mat(locations.t(), vValues, iNq*iNr*iNu,iNp*iNt*iNs);
        Nqru_pts = 1;
        return Vqru_pts;
    }
    else{
        return Vqru_pts;
    }
}

void flexmat6::update_as_qru_spt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vu*iNq*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vt*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qru_spt(){
    if(Nqru_spt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vu*iNr*iNq;
        locations.col(1) = vs + vp*iNs + vt*iNp*iNs;
        Vqru_spt = sp_mat(locations.t(), vValues, iNq*iNr*iNu,iNs*iNp*iNt);
        Nqru_spt = 1;
        return Vqru_spt;
    }
    else{
        return Vqru_spt;
    }
}

void flexmat6::update_as_qru_stp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vu*iNq*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vp*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qru_stp(){
    if(Nqru_stp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vu*iNr*iNq;
        locations.col(1) = vs + vt*iNs + vp*iNt*iNs;
        Vqru_stp = sp_mat(locations.t(), vValues, iNq*iNr*iNu,iNs*iNt*iNp);
        Nqru_stp = 1;
        return Vqru_stp;
    }
    else{
        return Vqru_stp;
    }
}

void flexmat6::update_as_qru_tps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vu*iNq*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vs*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qru_tps(){
    if(Nqru_tps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vu*iNr*iNq;
        locations.col(1) = vt + vp*iNt + vs*iNp*iNt;
        Vqru_tps = sp_mat(locations.t(), vValues, iNq*iNr*iNu,iNt*iNp*iNs);
        Nqru_tps = 1;
        return Vqru_tps;
    }
    else{
        return Vqru_tps;
    }
}

void flexmat6::update_as_qru_tsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNq - vu*iNq*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vp*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qru_tsp(){
    if(Nqru_tsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vr*iNq + vu*iNr*iNq;
        locations.col(1) = vt + vs*iNt + vp*iNs*iNt;
        Vqru_tsp = sp_mat(locations.t(), vValues, iNq*iNr*iNu,iNt*iNs*iNp);
        Nqru_tsp = 1;
        return Vqru_tsp;
    }
    else{
        return Vqru_tsp;
    }
}

void flexmat6::update_as_qsp_rtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vp*iNq*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vu*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsp_rtu(){
    if(Nqsp_rtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vp*iNs*iNq;
        locations.col(1) = vr + vt*iNr + vu*iNt*iNr;
        Vqsp_rtu = sp_mat(locations.t(), vValues, iNq*iNs*iNp,iNr*iNt*iNu);
        Nqsp_rtu = 1;
        return Vqsp_rtu;
    }
    else{
        return Vqsp_rtu;
    }
}

void flexmat6::update_as_qsp_rut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vp*iNq*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vt*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsp_rut(){
    if(Nqsp_rut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vp*iNs*iNq;
        locations.col(1) = vr + vu*iNr + vt*iNu*iNr;
        Vqsp_rut = sp_mat(locations.t(), vValues, iNq*iNs*iNp,iNr*iNu*iNt);
        Nqsp_rut = 1;
        return Vqsp_rut;
    }
    else{
        return Vqsp_rut;
    }
}

void flexmat6::update_as_qsp_tru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vp*iNq*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vu*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsp_tru(){
    if(Nqsp_tru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vp*iNs*iNq;
        locations.col(1) = vt + vr*iNt + vu*iNr*iNt;
        Vqsp_tru = sp_mat(locations.t(), vValues, iNq*iNs*iNp,iNt*iNr*iNu);
        Nqsp_tru = 1;
        return Vqsp_tru;
    }
    else{
        return Vqsp_tru;
    }
}

void flexmat6::update_as_qsp_tur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vp*iNq*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vr*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsp_tur(){
    if(Nqsp_tur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vp*iNs*iNq;
        locations.col(1) = vt + vu*iNt + vr*iNu*iNt;
        Vqsp_tur = sp_mat(locations.t(), vValues, iNq*iNs*iNp,iNt*iNu*iNr);
        Nqsp_tur = 1;
        return Vqsp_tur;
    }
    else{
        return Vqsp_tur;
    }
}

void flexmat6::update_as_qsp_urt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vp*iNq*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vt*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsp_urt(){
    if(Nqsp_urt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vp*iNs*iNq;
        locations.col(1) = vu + vr*iNu + vt*iNr*iNu;
        Vqsp_urt = sp_mat(locations.t(), vValues, iNq*iNs*iNp,iNu*iNr*iNt);
        Nqsp_urt = 1;
        return Vqsp_urt;
    }
    else{
        return Vqsp_urt;
    }
}

void flexmat6::update_as_qsp_utr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vp*iNq*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vr*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsp_utr(){
    if(Nqsp_utr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vp*iNs*iNq;
        locations.col(1) = vu + vt*iNu + vr*iNt*iNu;
        Vqsp_utr = sp_mat(locations.t(), vValues, iNq*iNs*iNp,iNu*iNt*iNr);
        Nqsp_utr = 1;
        return Vqsp_utr;
    }
    else{
        return Vqsp_utr;
    }
}

void flexmat6::update_as_qsr_ptu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vr*iNq*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vu*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsr_ptu(){
    if(Nqsr_ptu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vr*iNs*iNq;
        locations.col(1) = vp + vt*iNp + vu*iNt*iNp;
        Vqsr_ptu = sp_mat(locations.t(), vValues, iNq*iNs*iNr,iNp*iNt*iNu);
        Nqsr_ptu = 1;
        return Vqsr_ptu;
    }
    else{
        return Vqsr_ptu;
    }
}

void flexmat6::update_as_qsr_put(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vr*iNq*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vt*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsr_put(){
    if(Nqsr_put == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vr*iNs*iNq;
        locations.col(1) = vp + vu*iNp + vt*iNu*iNp;
        Vqsr_put = sp_mat(locations.t(), vValues, iNq*iNs*iNr,iNp*iNu*iNt);
        Nqsr_put = 1;
        return Vqsr_put;
    }
    else{
        return Vqsr_put;
    }
}

void flexmat6::update_as_qsr_tpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vr*iNq*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vu*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsr_tpu(){
    if(Nqsr_tpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vr*iNs*iNq;
        locations.col(1) = vt + vp*iNt + vu*iNp*iNt;
        Vqsr_tpu = sp_mat(locations.t(), vValues, iNq*iNs*iNr,iNt*iNp*iNu);
        Nqsr_tpu = 1;
        return Vqsr_tpu;
    }
    else{
        return Vqsr_tpu;
    }
}

void flexmat6::update_as_qsr_tup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vr*iNq*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vp*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsr_tup(){
    if(Nqsr_tup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vr*iNs*iNq;
        locations.col(1) = vt + vu*iNt + vp*iNu*iNt;
        Vqsr_tup = sp_mat(locations.t(), vValues, iNq*iNs*iNr,iNt*iNu*iNp);
        Nqsr_tup = 1;
        return Vqsr_tup;
    }
    else{
        return Vqsr_tup;
    }
}

void flexmat6::update_as_qsr_upt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vr*iNq*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vt*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsr_upt(){
    if(Nqsr_upt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vr*iNs*iNq;
        locations.col(1) = vu + vp*iNu + vt*iNp*iNu;
        Vqsr_upt = sp_mat(locations.t(), vValues, iNq*iNs*iNr,iNu*iNp*iNt);
        Nqsr_upt = 1;
        return Vqsr_upt;
    }
    else{
        return Vqsr_upt;
    }
}

void flexmat6::update_as_qsr_utp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vr*iNq*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vp*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsr_utp(){
    if(Nqsr_utp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vr*iNs*iNq;
        locations.col(1) = vu + vt*iNu + vp*iNt*iNu;
        Vqsr_utp = sp_mat(locations.t(), vValues, iNq*iNs*iNr,iNu*iNt*iNp);
        Nqsr_utp = 1;
        return Vqsr_utp;
    }
    else{
        return Vqsr_utp;
    }
}

void flexmat6::update_as_qst_pru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vt*iNq*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vu*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qst_pru(){
    if(Nqst_pru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vt*iNs*iNq;
        locations.col(1) = vp + vr*iNp + vu*iNr*iNp;
        Vqst_pru = sp_mat(locations.t(), vValues, iNq*iNs*iNt,iNp*iNr*iNu);
        Nqst_pru = 1;
        return Vqst_pru;
    }
    else{
        return Vqst_pru;
    }
}

void flexmat6::update_as_qst_pur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vt*iNq*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vr*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qst_pur(){
    if(Nqst_pur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vt*iNs*iNq;
        locations.col(1) = vp + vu*iNp + vr*iNu*iNp;
        Vqst_pur = sp_mat(locations.t(), vValues, iNq*iNs*iNt,iNp*iNu*iNr);
        Nqst_pur = 1;
        return Vqst_pur;
    }
    else{
        return Vqst_pur;
    }
}

void flexmat6::update_as_qst_rpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vt*iNq*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vu*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qst_rpu(){
    if(Nqst_rpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vt*iNs*iNq;
        locations.col(1) = vr + vp*iNr + vu*iNp*iNr;
        Vqst_rpu = sp_mat(locations.t(), vValues, iNq*iNs*iNt,iNr*iNp*iNu);
        Nqst_rpu = 1;
        return Vqst_rpu;
    }
    else{
        return Vqst_rpu;
    }
}

void flexmat6::update_as_qst_rup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vt*iNq*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vp*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qst_rup(){
    if(Nqst_rup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vt*iNs*iNq;
        locations.col(1) = vr + vu*iNr + vp*iNu*iNr;
        Vqst_rup = sp_mat(locations.t(), vValues, iNq*iNs*iNt,iNr*iNu*iNp);
        Nqst_rup = 1;
        return Vqst_rup;
    }
    else{
        return Vqst_rup;
    }
}

void flexmat6::update_as_qst_upr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vt*iNq*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vr*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qst_upr(){
    if(Nqst_upr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vt*iNs*iNq;
        locations.col(1) = vu + vp*iNu + vr*iNp*iNu;
        Vqst_upr = sp_mat(locations.t(), vValues, iNq*iNs*iNt,iNu*iNp*iNr);
        Nqst_upr = 1;
        return Vqst_upr;
    }
    else{
        return Vqst_upr;
    }
}

void flexmat6::update_as_qst_urp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vt*iNq*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vp*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qst_urp(){
    if(Nqst_urp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vt*iNs*iNq;
        locations.col(1) = vu + vr*iNu + vp*iNr*iNu;
        Vqst_urp = sp_mat(locations.t(), vValues, iNq*iNs*iNt,iNu*iNr*iNp);
        Nqst_urp = 1;
        return Vqst_urp;
    }
    else{
        return Vqst_urp;
    }
}

void flexmat6::update_as_qsu_prt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vu*iNq*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vt*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsu_prt(){
    if(Nqsu_prt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vu*iNs*iNq;
        locations.col(1) = vp + vr*iNp + vt*iNr*iNp;
        Vqsu_prt = sp_mat(locations.t(), vValues, iNq*iNs*iNu,iNp*iNr*iNt);
        Nqsu_prt = 1;
        return Vqsu_prt;
    }
    else{
        return Vqsu_prt;
    }
}

void flexmat6::update_as_qsu_ptr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vu*iNq*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vr*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsu_ptr(){
    if(Nqsu_ptr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vu*iNs*iNq;
        locations.col(1) = vp + vt*iNp + vr*iNt*iNp;
        Vqsu_ptr = sp_mat(locations.t(), vValues, iNq*iNs*iNu,iNp*iNt*iNr);
        Nqsu_ptr = 1;
        return Vqsu_ptr;
    }
    else{
        return Vqsu_ptr;
    }
}

void flexmat6::update_as_qsu_rpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vu*iNq*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vt*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsu_rpt(){
    if(Nqsu_rpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vu*iNs*iNq;
        locations.col(1) = vr + vp*iNr + vt*iNp*iNr;
        Vqsu_rpt = sp_mat(locations.t(), vValues, iNq*iNs*iNu,iNr*iNp*iNt);
        Nqsu_rpt = 1;
        return Vqsu_rpt;
    }
    else{
        return Vqsu_rpt;
    }
}

void flexmat6::update_as_qsu_rtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vu*iNq*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vp*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsu_rtp(){
    if(Nqsu_rtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vu*iNs*iNq;
        locations.col(1) = vr + vt*iNr + vp*iNt*iNr;
        Vqsu_rtp = sp_mat(locations.t(), vValues, iNq*iNs*iNu,iNr*iNt*iNp);
        Nqsu_rtp = 1;
        return Vqsu_rtp;
    }
    else{
        return Vqsu_rtp;
    }
}

void flexmat6::update_as_qsu_tpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vu*iNq*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vr*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsu_tpr(){
    if(Nqsu_tpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vu*iNs*iNq;
        locations.col(1) = vt + vp*iNt + vr*iNp*iNt;
        Vqsu_tpr = sp_mat(locations.t(), vValues, iNq*iNs*iNu,iNt*iNp*iNr);
        Nqsu_tpr = 1;
        return Vqsu_tpr;
    }
    else{
        return Vqsu_tpr;
    }
}

void flexmat6::update_as_qsu_trp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNq - vu*iNq*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vp*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qsu_trp(){
    if(Nqsu_trp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vs*iNq + vu*iNs*iNq;
        locations.col(1) = vt + vr*iNt + vp*iNr*iNt;
        Vqsu_trp = sp_mat(locations.t(), vValues, iNq*iNs*iNu,iNt*iNr*iNp);
        Nqsu_trp = 1;
        return Vqsu_trp;
    }
    else{
        return Vqsu_trp;
    }
}

void flexmat6::update_as_qtp_rsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vp*iNq*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vu*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtp_rsu(){
    if(Nqtp_rsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vp*iNt*iNq;
        locations.col(1) = vr + vs*iNr + vu*iNs*iNr;
        Vqtp_rsu = sp_mat(locations.t(), vValues, iNq*iNt*iNp,iNr*iNs*iNu);
        Nqtp_rsu = 1;
        return Vqtp_rsu;
    }
    else{
        return Vqtp_rsu;
    }
}

void flexmat6::update_as_qtp_rus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vp*iNq*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vs*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtp_rus(){
    if(Nqtp_rus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vp*iNt*iNq;
        locations.col(1) = vr + vu*iNr + vs*iNu*iNr;
        Vqtp_rus = sp_mat(locations.t(), vValues, iNq*iNt*iNp,iNr*iNu*iNs);
        Nqtp_rus = 1;
        return Vqtp_rus;
    }
    else{
        return Vqtp_rus;
    }
}

void flexmat6::update_as_qtp_sru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vp*iNq*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vu*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtp_sru(){
    if(Nqtp_sru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vp*iNt*iNq;
        locations.col(1) = vs + vr*iNs + vu*iNr*iNs;
        Vqtp_sru = sp_mat(locations.t(), vValues, iNq*iNt*iNp,iNs*iNr*iNu);
        Nqtp_sru = 1;
        return Vqtp_sru;
    }
    else{
        return Vqtp_sru;
    }
}

void flexmat6::update_as_qtp_sur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vp*iNq*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vr*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtp_sur(){
    if(Nqtp_sur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vp*iNt*iNq;
        locations.col(1) = vs + vu*iNs + vr*iNu*iNs;
        Vqtp_sur = sp_mat(locations.t(), vValues, iNq*iNt*iNp,iNs*iNu*iNr);
        Nqtp_sur = 1;
        return Vqtp_sur;
    }
    else{
        return Vqtp_sur;
    }
}

void flexmat6::update_as_qtp_urs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vp*iNq*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vs*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtp_urs(){
    if(Nqtp_urs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vp*iNt*iNq;
        locations.col(1) = vu + vr*iNu + vs*iNr*iNu;
        Vqtp_urs = sp_mat(locations.t(), vValues, iNq*iNt*iNp,iNu*iNr*iNs);
        Nqtp_urs = 1;
        return Vqtp_urs;
    }
    else{
        return Vqtp_urs;
    }
}

void flexmat6::update_as_qtp_usr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vp*iNq*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vr*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtp_usr(){
    if(Nqtp_usr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vp*iNt*iNq;
        locations.col(1) = vu + vs*iNu + vr*iNs*iNu;
        Vqtp_usr = sp_mat(locations.t(), vValues, iNq*iNt*iNp,iNu*iNs*iNr);
        Nqtp_usr = 1;
        return Vqtp_usr;
    }
    else{
        return Vqtp_usr;
    }
}

void flexmat6::update_as_qtr_psu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vr*iNq*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vu*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtr_psu(){
    if(Nqtr_psu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vr*iNt*iNq;
        locations.col(1) = vp + vs*iNp + vu*iNs*iNp;
        Vqtr_psu = sp_mat(locations.t(), vValues, iNq*iNt*iNr,iNp*iNs*iNu);
        Nqtr_psu = 1;
        return Vqtr_psu;
    }
    else{
        return Vqtr_psu;
    }
}

void flexmat6::update_as_qtr_pus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vr*iNq*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vs*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtr_pus(){
    if(Nqtr_pus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vr*iNt*iNq;
        locations.col(1) = vp + vu*iNp + vs*iNu*iNp;
        Vqtr_pus = sp_mat(locations.t(), vValues, iNq*iNt*iNr,iNp*iNu*iNs);
        Nqtr_pus = 1;
        return Vqtr_pus;
    }
    else{
        return Vqtr_pus;
    }
}

void flexmat6::update_as_qtr_spu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vr*iNq*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vu*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtr_spu(){
    if(Nqtr_spu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vr*iNt*iNq;
        locations.col(1) = vs + vp*iNs + vu*iNp*iNs;
        Vqtr_spu = sp_mat(locations.t(), vValues, iNq*iNt*iNr,iNs*iNp*iNu);
        Nqtr_spu = 1;
        return Vqtr_spu;
    }
    else{
        return Vqtr_spu;
    }
}

void flexmat6::update_as_qtr_sup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vr*iNq*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vp*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtr_sup(){
    if(Nqtr_sup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vr*iNt*iNq;
        locations.col(1) = vs + vu*iNs + vp*iNu*iNs;
        Vqtr_sup = sp_mat(locations.t(), vValues, iNq*iNt*iNr,iNs*iNu*iNp);
        Nqtr_sup = 1;
        return Vqtr_sup;
    }
    else{
        return Vqtr_sup;
    }
}

void flexmat6::update_as_qtr_ups(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vr*iNq*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vs*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtr_ups(){
    if(Nqtr_ups == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vr*iNt*iNq;
        locations.col(1) = vu + vp*iNu + vs*iNp*iNu;
        Vqtr_ups = sp_mat(locations.t(), vValues, iNq*iNt*iNr,iNu*iNp*iNs);
        Nqtr_ups = 1;
        return Vqtr_ups;
    }
    else{
        return Vqtr_ups;
    }
}

void flexmat6::update_as_qtr_usp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vr*iNq*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vp*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtr_usp(){
    if(Nqtr_usp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vr*iNt*iNq;
        locations.col(1) = vu + vs*iNu + vp*iNs*iNu;
        Vqtr_usp = sp_mat(locations.t(), vValues, iNq*iNt*iNr,iNu*iNs*iNp);
        Nqtr_usp = 1;
        return Vqtr_usp;
    }
    else{
        return Vqtr_usp;
    }
}

void flexmat6::update_as_qts_pru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vs*iNq*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vu*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qts_pru(){
    if(Nqts_pru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vs*iNt*iNq;
        locations.col(1) = vp + vr*iNp + vu*iNr*iNp;
        Vqts_pru = sp_mat(locations.t(), vValues, iNq*iNt*iNs,iNp*iNr*iNu);
        Nqts_pru = 1;
        return Vqts_pru;
    }
    else{
        return Vqts_pru;
    }
}

void flexmat6::update_as_qts_pur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vs*iNq*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vr*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qts_pur(){
    if(Nqts_pur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vs*iNt*iNq;
        locations.col(1) = vp + vu*iNp + vr*iNu*iNp;
        Vqts_pur = sp_mat(locations.t(), vValues, iNq*iNt*iNs,iNp*iNu*iNr);
        Nqts_pur = 1;
        return Vqts_pur;
    }
    else{
        return Vqts_pur;
    }
}

void flexmat6::update_as_qts_rpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vs*iNq*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vu*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qts_rpu(){
    if(Nqts_rpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vs*iNt*iNq;
        locations.col(1) = vr + vp*iNr + vu*iNp*iNr;
        Vqts_rpu = sp_mat(locations.t(), vValues, iNq*iNt*iNs,iNr*iNp*iNu);
        Nqts_rpu = 1;
        return Vqts_rpu;
    }
    else{
        return Vqts_rpu;
    }
}

void flexmat6::update_as_qts_rup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vs*iNq*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vp*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qts_rup(){
    if(Nqts_rup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vs*iNt*iNq;
        locations.col(1) = vr + vu*iNr + vp*iNu*iNr;
        Vqts_rup = sp_mat(locations.t(), vValues, iNq*iNt*iNs,iNr*iNu*iNp);
        Nqts_rup = 1;
        return Vqts_rup;
    }
    else{
        return Vqts_rup;
    }
}

void flexmat6::update_as_qts_upr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vs*iNq*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vr*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qts_upr(){
    if(Nqts_upr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vs*iNt*iNq;
        locations.col(1) = vu + vp*iNu + vr*iNp*iNu;
        Vqts_upr = sp_mat(locations.t(), vValues, iNq*iNt*iNs,iNu*iNp*iNr);
        Nqts_upr = 1;
        return Vqts_upr;
    }
    else{
        return Vqts_upr;
    }
}

void flexmat6::update_as_qts_urp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vs*iNq*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vp*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qts_urp(){
    if(Nqts_urp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vs*iNt*iNq;
        locations.col(1) = vu + vr*iNu + vp*iNr*iNu;
        Vqts_urp = sp_mat(locations.t(), vValues, iNq*iNt*iNs,iNu*iNr*iNp);
        Nqts_urp = 1;
        return Vqts_urp;
    }
    else{
        return Vqts_urp;
    }
}

void flexmat6::update_as_qtu_prs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vu*iNq*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vs*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtu_prs(){
    if(Nqtu_prs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vu*iNt*iNq;
        locations.col(1) = vp + vr*iNp + vs*iNr*iNp;
        Vqtu_prs = sp_mat(locations.t(), vValues, iNq*iNt*iNu,iNp*iNr*iNs);
        Nqtu_prs = 1;
        return Vqtu_prs;
    }
    else{
        return Vqtu_prs;
    }
}

void flexmat6::update_as_qtu_psr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vu*iNq*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vr*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtu_psr(){
    if(Nqtu_psr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vu*iNt*iNq;
        locations.col(1) = vp + vs*iNp + vr*iNs*iNp;
        Vqtu_psr = sp_mat(locations.t(), vValues, iNq*iNt*iNu,iNp*iNs*iNr);
        Nqtu_psr = 1;
        return Vqtu_psr;
    }
    else{
        return Vqtu_psr;
    }
}

void flexmat6::update_as_qtu_rps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vu*iNq*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vs*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtu_rps(){
    if(Nqtu_rps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vu*iNt*iNq;
        locations.col(1) = vr + vp*iNr + vs*iNp*iNr;
        Vqtu_rps = sp_mat(locations.t(), vValues, iNq*iNt*iNu,iNr*iNp*iNs);
        Nqtu_rps = 1;
        return Vqtu_rps;
    }
    else{
        return Vqtu_rps;
    }
}

void flexmat6::update_as_qtu_rsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vu*iNq*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vp*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtu_rsp(){
    if(Nqtu_rsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vu*iNt*iNq;
        locations.col(1) = vr + vs*iNr + vp*iNs*iNr;
        Vqtu_rsp = sp_mat(locations.t(), vValues, iNq*iNt*iNu,iNr*iNs*iNp);
        Nqtu_rsp = 1;
        return Vqtu_rsp;
    }
    else{
        return Vqtu_rsp;
    }
}

void flexmat6::update_as_qtu_spr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vu*iNq*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vr*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtu_spr(){
    if(Nqtu_spr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vu*iNt*iNq;
        locations.col(1) = vs + vp*iNs + vr*iNp*iNs;
        Vqtu_spr = sp_mat(locations.t(), vValues, iNq*iNt*iNu,iNs*iNp*iNr);
        Nqtu_spr = 1;
        return Vqtu_spr;
    }
    else{
        return Vqtu_spr;
    }
}

void flexmat6::update_as_qtu_srp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNq - vu*iNq*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vp*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qtu_srp(){
    if(Nqtu_srp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vt*iNq + vu*iNt*iNq;
        locations.col(1) = vs + vr*iNs + vp*iNr*iNs;
        Vqtu_srp = sp_mat(locations.t(), vValues, iNq*iNt*iNu,iNs*iNr*iNp);
        Nqtu_srp = 1;
        return Vqtu_srp;
    }
    else{
        return Vqtu_srp;
    }
}

void flexmat6::update_as_qup_rst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vp*iNq*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vt*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qup_rst(){
    if(Nqup_rst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vp*iNu*iNq;
        locations.col(1) = vr + vs*iNr + vt*iNs*iNr;
        Vqup_rst = sp_mat(locations.t(), vValues, iNq*iNu*iNp,iNr*iNs*iNt);
        Nqup_rst = 1;
        return Vqup_rst;
    }
    else{
        return Vqup_rst;
    }
}

void flexmat6::update_as_qup_rts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vp*iNq*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vs*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qup_rts(){
    if(Nqup_rts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vp*iNu*iNq;
        locations.col(1) = vr + vt*iNr + vs*iNt*iNr;
        Vqup_rts = sp_mat(locations.t(), vValues, iNq*iNu*iNp,iNr*iNt*iNs);
        Nqup_rts = 1;
        return Vqup_rts;
    }
    else{
        return Vqup_rts;
    }
}

void flexmat6::update_as_qup_srt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vp*iNq*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vt*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qup_srt(){
    if(Nqup_srt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vp*iNu*iNq;
        locations.col(1) = vs + vr*iNs + vt*iNr*iNs;
        Vqup_srt = sp_mat(locations.t(), vValues, iNq*iNu*iNp,iNs*iNr*iNt);
        Nqup_srt = 1;
        return Vqup_srt;
    }
    else{
        return Vqup_srt;
    }
}

void flexmat6::update_as_qup_str(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vp*iNq*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vr*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qup_str(){
    if(Nqup_str == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vp*iNu*iNq;
        locations.col(1) = vs + vt*iNs + vr*iNt*iNs;
        Vqup_str = sp_mat(locations.t(), vValues, iNq*iNu*iNp,iNs*iNt*iNr);
        Nqup_str = 1;
        return Vqup_str;
    }
    else{
        return Vqup_str;
    }
}

void flexmat6::update_as_qup_trs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vp*iNq*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vs*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qup_trs(){
    if(Nqup_trs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vp*iNu*iNq;
        locations.col(1) = vt + vr*iNt + vs*iNr*iNt;
        Vqup_trs = sp_mat(locations.t(), vValues, iNq*iNu*iNp,iNt*iNr*iNs);
        Nqup_trs = 1;
        return Vqup_trs;
    }
    else{
        return Vqup_trs;
    }
}

void flexmat6::update_as_qup_tsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vp*iNq*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vr*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qup_tsr(){
    if(Nqup_tsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vp*iNu*iNq;
        locations.col(1) = vt + vs*iNt + vr*iNs*iNt;
        Vqup_tsr = sp_mat(locations.t(), vValues, iNq*iNu*iNp,iNt*iNs*iNr);
        Nqup_tsr = 1;
        return Vqup_tsr;
    }
    else{
        return Vqup_tsr;
    }
}

void flexmat6::update_as_qur_pst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vr*iNq*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vt*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qur_pst(){
    if(Nqur_pst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vr*iNu*iNq;
        locations.col(1) = vp + vs*iNp + vt*iNs*iNp;
        Vqur_pst = sp_mat(locations.t(), vValues, iNq*iNu*iNr,iNp*iNs*iNt);
        Nqur_pst = 1;
        return Vqur_pst;
    }
    else{
        return Vqur_pst;
    }
}

void flexmat6::update_as_qur_pts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vr*iNq*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vs*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qur_pts(){
    if(Nqur_pts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vr*iNu*iNq;
        locations.col(1) = vp + vt*iNp + vs*iNt*iNp;
        Vqur_pts = sp_mat(locations.t(), vValues, iNq*iNu*iNr,iNp*iNt*iNs);
        Nqur_pts = 1;
        return Vqur_pts;
    }
    else{
        return Vqur_pts;
    }
}

void flexmat6::update_as_qur_spt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vr*iNq*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vt*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qur_spt(){
    if(Nqur_spt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vr*iNu*iNq;
        locations.col(1) = vs + vp*iNs + vt*iNp*iNs;
        Vqur_spt = sp_mat(locations.t(), vValues, iNq*iNu*iNr,iNs*iNp*iNt);
        Nqur_spt = 1;
        return Vqur_spt;
    }
    else{
        return Vqur_spt;
    }
}

void flexmat6::update_as_qur_stp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vr*iNq*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vp*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qur_stp(){
    if(Nqur_stp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vr*iNu*iNq;
        locations.col(1) = vs + vt*iNs + vp*iNt*iNs;
        Vqur_stp = sp_mat(locations.t(), vValues, iNq*iNu*iNr,iNs*iNt*iNp);
        Nqur_stp = 1;
        return Vqur_stp;
    }
    else{
        return Vqur_stp;
    }
}

void flexmat6::update_as_qur_tps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vr*iNq*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vs*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qur_tps(){
    if(Nqur_tps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vr*iNu*iNq;
        locations.col(1) = vt + vp*iNt + vs*iNp*iNt;
        Vqur_tps = sp_mat(locations.t(), vValues, iNq*iNu*iNr,iNt*iNp*iNs);
        Nqur_tps = 1;
        return Vqur_tps;
    }
    else{
        return Vqur_tps;
    }
}

void flexmat6::update_as_qur_tsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vr*iNq*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vp*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qur_tsp(){
    if(Nqur_tsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vr*iNu*iNq;
        locations.col(1) = vt + vs*iNt + vp*iNs*iNt;
        Vqur_tsp = sp_mat(locations.t(), vValues, iNq*iNu*iNr,iNt*iNs*iNp);
        Nqur_tsp = 1;
        return Vqur_tsp;
    }
    else{
        return Vqur_tsp;
    }
}

void flexmat6::update_as_qus_prt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vs*iNq*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vt*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qus_prt(){
    if(Nqus_prt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vs*iNu*iNq;
        locations.col(1) = vp + vr*iNp + vt*iNr*iNp;
        Vqus_prt = sp_mat(locations.t(), vValues, iNq*iNu*iNs,iNp*iNr*iNt);
        Nqus_prt = 1;
        return Vqus_prt;
    }
    else{
        return Vqus_prt;
    }
}

void flexmat6::update_as_qus_ptr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vs*iNq*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vr*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qus_ptr(){
    if(Nqus_ptr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vs*iNu*iNq;
        locations.col(1) = vp + vt*iNp + vr*iNt*iNp;
        Vqus_ptr = sp_mat(locations.t(), vValues, iNq*iNu*iNs,iNp*iNt*iNr);
        Nqus_ptr = 1;
        return Vqus_ptr;
    }
    else{
        return Vqus_ptr;
    }
}

void flexmat6::update_as_qus_rpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vs*iNq*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vt*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qus_rpt(){
    if(Nqus_rpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vs*iNu*iNq;
        locations.col(1) = vr + vp*iNr + vt*iNp*iNr;
        Vqus_rpt = sp_mat(locations.t(), vValues, iNq*iNu*iNs,iNr*iNp*iNt);
        Nqus_rpt = 1;
        return Vqus_rpt;
    }
    else{
        return Vqus_rpt;
    }
}

void flexmat6::update_as_qus_rtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vs*iNq*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vp*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qus_rtp(){
    if(Nqus_rtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vs*iNu*iNq;
        locations.col(1) = vr + vt*iNr + vp*iNt*iNr;
        Vqus_rtp = sp_mat(locations.t(), vValues, iNq*iNu*iNs,iNr*iNt*iNp);
        Nqus_rtp = 1;
        return Vqus_rtp;
    }
    else{
        return Vqus_rtp;
    }
}

void flexmat6::update_as_qus_tpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vs*iNq*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vr*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qus_tpr(){
    if(Nqus_tpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vs*iNu*iNq;
        locations.col(1) = vt + vp*iNt + vr*iNp*iNt;
        Vqus_tpr = sp_mat(locations.t(), vValues, iNq*iNu*iNs,iNt*iNp*iNr);
        Nqus_tpr = 1;
        return Vqus_tpr;
    }
    else{
        return Vqus_tpr;
    }
}

void flexmat6::update_as_qus_trp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vs*iNq*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vp*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qus_trp(){
    if(Nqus_trp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vs*iNu*iNq;
        locations.col(1) = vt + vr*iNt + vp*iNr*iNt;
        Vqus_trp = sp_mat(locations.t(), vValues, iNq*iNu*iNs,iNt*iNr*iNp);
        Nqus_trp = 1;
        return Vqus_trp;
    }
    else{
        return Vqus_trp;
    }
}

void flexmat6::update_as_qut_prs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vt*iNq*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vs*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qut_prs(){
    if(Nqut_prs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vt*iNu*iNq;
        locations.col(1) = vp + vr*iNp + vs*iNr*iNp;
        Vqut_prs = sp_mat(locations.t(), vValues, iNq*iNu*iNt,iNp*iNr*iNs);
        Nqut_prs = 1;
        return Vqut_prs;
    }
    else{
        return Vqut_prs;
    }
}

void flexmat6::update_as_qut_psr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vt*iNq*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vr*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qut_psr(){
    if(Nqut_psr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vt*iNu*iNq;
        locations.col(1) = vp + vs*iNp + vr*iNs*iNp;
        Vqut_psr = sp_mat(locations.t(), vValues, iNq*iNu*iNt,iNp*iNs*iNr);
        Nqut_psr = 1;
        return Vqut_psr;
    }
    else{
        return Vqut_psr;
    }
}

void flexmat6::update_as_qut_rps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vt*iNq*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vs*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qut_rps(){
    if(Nqut_rps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vt*iNu*iNq;
        locations.col(1) = vr + vp*iNr + vs*iNp*iNr;
        Vqut_rps = sp_mat(locations.t(), vValues, iNq*iNu*iNt,iNr*iNp*iNs);
        Nqut_rps = 1;
        return Vqut_rps;
    }
    else{
        return Vqut_rps;
    }
}

void flexmat6::update_as_qut_rsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vt*iNq*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vp*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qut_rsp(){
    if(Nqut_rsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vt*iNu*iNq;
        locations.col(1) = vr + vs*iNr + vp*iNs*iNr;
        Vqut_rsp = sp_mat(locations.t(), vValues, iNq*iNu*iNt,iNr*iNs*iNp);
        Nqut_rsp = 1;
        return Vqut_rsp;
    }
    else{
        return Vqut_rsp;
    }
}

void flexmat6::update_as_qut_spr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vt*iNq*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vr*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qut_spr(){
    if(Nqut_spr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vt*iNu*iNq;
        locations.col(1) = vs + vp*iNs + vr*iNp*iNs;
        Vqut_spr = sp_mat(locations.t(), vValues, iNq*iNu*iNt,iNs*iNp*iNr);
        Nqut_spr = 1;
        return Vqut_spr;
    }
    else{
        return Vqut_spr;
    }
}

void flexmat6::update_as_qut_srp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNq - vt*iNq*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vp*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::qut_srp(){
    if(Nqut_srp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vq + vu*iNq + vt*iNu*iNq;
        locations.col(1) = vs + vr*iNs + vp*iNr*iNs;
        Vqut_srp = sp_mat(locations.t(), vValues, iNq*iNu*iNt,iNs*iNr*iNp);
        Nqut_srp = 1;
        return Vqut_srp;
    }
    else{
        return Vqut_srp;
    }
}

void flexmat6::update_as_rpq_stu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vq*iNr*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vu*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpq_stu(){
    if(Nrpq_stu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vq*iNp*iNr;
        locations.col(1) = vs + vt*iNs + vu*iNt*iNs;
        Vrpq_stu = sp_mat(locations.t(), vValues, iNr*iNp*iNq,iNs*iNt*iNu);
        Nrpq_stu = 1;
        return Vrpq_stu;
    }
    else{
        return Vrpq_stu;
    }
}

void flexmat6::update_as_rpq_sut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vq*iNr*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vt*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpq_sut(){
    if(Nrpq_sut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vq*iNp*iNr;
        locations.col(1) = vs + vu*iNs + vt*iNu*iNs;
        Vrpq_sut = sp_mat(locations.t(), vValues, iNr*iNp*iNq,iNs*iNu*iNt);
        Nrpq_sut = 1;
        return Vrpq_sut;
    }
    else{
        return Vrpq_sut;
    }
}

void flexmat6::update_as_rpq_tsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vq*iNr*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vu*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpq_tsu(){
    if(Nrpq_tsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vq*iNp*iNr;
        locations.col(1) = vt + vs*iNt + vu*iNs*iNt;
        Vrpq_tsu = sp_mat(locations.t(), vValues, iNr*iNp*iNq,iNt*iNs*iNu);
        Nrpq_tsu = 1;
        return Vrpq_tsu;
    }
    else{
        return Vrpq_tsu;
    }
}

void flexmat6::update_as_rpq_tus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vq*iNr*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vs*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpq_tus(){
    if(Nrpq_tus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vq*iNp*iNr;
        locations.col(1) = vt + vu*iNt + vs*iNu*iNt;
        Vrpq_tus = sp_mat(locations.t(), vValues, iNr*iNp*iNq,iNt*iNu*iNs);
        Nrpq_tus = 1;
        return Vrpq_tus;
    }
    else{
        return Vrpq_tus;
    }
}

void flexmat6::update_as_rpq_ust(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vq*iNr*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vt*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpq_ust(){
    if(Nrpq_ust == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vq*iNp*iNr;
        locations.col(1) = vu + vs*iNu + vt*iNs*iNu;
        Vrpq_ust = sp_mat(locations.t(), vValues, iNr*iNp*iNq,iNu*iNs*iNt);
        Nrpq_ust = 1;
        return Vrpq_ust;
    }
    else{
        return Vrpq_ust;
    }
}

void flexmat6::update_as_rpq_uts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vq*iNr*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vs*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpq_uts(){
    if(Nrpq_uts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vq*iNp*iNr;
        locations.col(1) = vu + vt*iNu + vs*iNt*iNu;
        Vrpq_uts = sp_mat(locations.t(), vValues, iNr*iNp*iNq,iNu*iNt*iNs);
        Nrpq_uts = 1;
        return Vrpq_uts;
    }
    else{
        return Vrpq_uts;
    }
}

void flexmat6::update_as_rps_qtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vs*iNr*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vu*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rps_qtu(){
    if(Nrps_qtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vs*iNp*iNr;
        locations.col(1) = vq + vt*iNq + vu*iNt*iNq;
        Vrps_qtu = sp_mat(locations.t(), vValues, iNr*iNp*iNs,iNq*iNt*iNu);
        Nrps_qtu = 1;
        return Vrps_qtu;
    }
    else{
        return Vrps_qtu;
    }
}

void flexmat6::update_as_rps_qut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vs*iNr*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vt*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rps_qut(){
    if(Nrps_qut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vs*iNp*iNr;
        locations.col(1) = vq + vu*iNq + vt*iNu*iNq;
        Vrps_qut = sp_mat(locations.t(), vValues, iNr*iNp*iNs,iNq*iNu*iNt);
        Nrps_qut = 1;
        return Vrps_qut;
    }
    else{
        return Vrps_qut;
    }
}

void flexmat6::update_as_rps_tqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vs*iNr*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vu*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rps_tqu(){
    if(Nrps_tqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vs*iNp*iNr;
        locations.col(1) = vt + vq*iNt + vu*iNq*iNt;
        Vrps_tqu = sp_mat(locations.t(), vValues, iNr*iNp*iNs,iNt*iNq*iNu);
        Nrps_tqu = 1;
        return Vrps_tqu;
    }
    else{
        return Vrps_tqu;
    }
}

void flexmat6::update_as_rps_tuq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vs*iNr*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vq*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rps_tuq(){
    if(Nrps_tuq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vs*iNp*iNr;
        locations.col(1) = vt + vu*iNt + vq*iNu*iNt;
        Vrps_tuq = sp_mat(locations.t(), vValues, iNr*iNp*iNs,iNt*iNu*iNq);
        Nrps_tuq = 1;
        return Vrps_tuq;
    }
    else{
        return Vrps_tuq;
    }
}

void flexmat6::update_as_rps_uqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vs*iNr*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vt*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rps_uqt(){
    if(Nrps_uqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vs*iNp*iNr;
        locations.col(1) = vu + vq*iNu + vt*iNq*iNu;
        Vrps_uqt = sp_mat(locations.t(), vValues, iNr*iNp*iNs,iNu*iNq*iNt);
        Nrps_uqt = 1;
        return Vrps_uqt;
    }
    else{
        return Vrps_uqt;
    }
}

void flexmat6::update_as_rps_utq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vs*iNr*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vq*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rps_utq(){
    if(Nrps_utq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vs*iNp*iNr;
        locations.col(1) = vu + vt*iNu + vq*iNt*iNu;
        Vrps_utq = sp_mat(locations.t(), vValues, iNr*iNp*iNs,iNu*iNt*iNq);
        Nrps_utq = 1;
        return Vrps_utq;
    }
    else{
        return Vrps_utq;
    }
}

void flexmat6::update_as_rpt_qsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vt*iNr*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vu*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpt_qsu(){
    if(Nrpt_qsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vt*iNp*iNr;
        locations.col(1) = vq + vs*iNq + vu*iNs*iNq;
        Vrpt_qsu = sp_mat(locations.t(), vValues, iNr*iNp*iNt,iNq*iNs*iNu);
        Nrpt_qsu = 1;
        return Vrpt_qsu;
    }
    else{
        return Vrpt_qsu;
    }
}

void flexmat6::update_as_rpt_qus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vt*iNr*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vs*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpt_qus(){
    if(Nrpt_qus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vt*iNp*iNr;
        locations.col(1) = vq + vu*iNq + vs*iNu*iNq;
        Vrpt_qus = sp_mat(locations.t(), vValues, iNr*iNp*iNt,iNq*iNu*iNs);
        Nrpt_qus = 1;
        return Vrpt_qus;
    }
    else{
        return Vrpt_qus;
    }
}

void flexmat6::update_as_rpt_squ(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vt*iNr*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vu*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpt_squ(){
    if(Nrpt_squ == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vt*iNp*iNr;
        locations.col(1) = vs + vq*iNs + vu*iNq*iNs;
        Vrpt_squ = sp_mat(locations.t(), vValues, iNr*iNp*iNt,iNs*iNq*iNu);
        Nrpt_squ = 1;
        return Vrpt_squ;
    }
    else{
        return Vrpt_squ;
    }
}

void flexmat6::update_as_rpt_suq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vt*iNr*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vq*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpt_suq(){
    if(Nrpt_suq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vt*iNp*iNr;
        locations.col(1) = vs + vu*iNs + vq*iNu*iNs;
        Vrpt_suq = sp_mat(locations.t(), vValues, iNr*iNp*iNt,iNs*iNu*iNq);
        Nrpt_suq = 1;
        return Vrpt_suq;
    }
    else{
        return Vrpt_suq;
    }
}

void flexmat6::update_as_rpt_uqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vt*iNr*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vs*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpt_uqs(){
    if(Nrpt_uqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vt*iNp*iNr;
        locations.col(1) = vu + vq*iNu + vs*iNq*iNu;
        Vrpt_uqs = sp_mat(locations.t(), vValues, iNr*iNp*iNt,iNu*iNq*iNs);
        Nrpt_uqs = 1;
        return Vrpt_uqs;
    }
    else{
        return Vrpt_uqs;
    }
}

void flexmat6::update_as_rpt_usq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vt*iNr*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vq*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpt_usq(){
    if(Nrpt_usq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vt*iNp*iNr;
        locations.col(1) = vu + vs*iNu + vq*iNs*iNu;
        Vrpt_usq = sp_mat(locations.t(), vValues, iNr*iNp*iNt,iNu*iNs*iNq);
        Nrpt_usq = 1;
        return Vrpt_usq;
    }
    else{
        return Vrpt_usq;
    }
}

void flexmat6::update_as_rpu_qst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vu*iNr*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vt*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpu_qst(){
    if(Nrpu_qst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vu*iNp*iNr;
        locations.col(1) = vq + vs*iNq + vt*iNs*iNq;
        Vrpu_qst = sp_mat(locations.t(), vValues, iNr*iNp*iNu,iNq*iNs*iNt);
        Nrpu_qst = 1;
        return Vrpu_qst;
    }
    else{
        return Vrpu_qst;
    }
}

void flexmat6::update_as_rpu_qts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vu*iNr*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vs*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpu_qts(){
    if(Nrpu_qts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vu*iNp*iNr;
        locations.col(1) = vq + vt*iNq + vs*iNt*iNq;
        Vrpu_qts = sp_mat(locations.t(), vValues, iNr*iNp*iNu,iNq*iNt*iNs);
        Nrpu_qts = 1;
        return Vrpu_qts;
    }
    else{
        return Vrpu_qts;
    }
}

void flexmat6::update_as_rpu_sqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vu*iNr*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vt*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpu_sqt(){
    if(Nrpu_sqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vu*iNp*iNr;
        locations.col(1) = vs + vq*iNs + vt*iNq*iNs;
        Vrpu_sqt = sp_mat(locations.t(), vValues, iNr*iNp*iNu,iNs*iNq*iNt);
        Nrpu_sqt = 1;
        return Vrpu_sqt;
    }
    else{
        return Vrpu_sqt;
    }
}

void flexmat6::update_as_rpu_stq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vu*iNr*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vq*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpu_stq(){
    if(Nrpu_stq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vu*iNp*iNr;
        locations.col(1) = vs + vt*iNs + vq*iNt*iNs;
        Vrpu_stq = sp_mat(locations.t(), vValues, iNr*iNp*iNu,iNs*iNt*iNq);
        Nrpu_stq = 1;
        return Vrpu_stq;
    }
    else{
        return Vrpu_stq;
    }
}

void flexmat6::update_as_rpu_tqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vu*iNr*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vs*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpu_tqs(){
    if(Nrpu_tqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vu*iNp*iNr;
        locations.col(1) = vt + vq*iNt + vs*iNq*iNt;
        Vrpu_tqs = sp_mat(locations.t(), vValues, iNr*iNp*iNu,iNt*iNq*iNs);
        Nrpu_tqs = 1;
        return Vrpu_tqs;
    }
    else{
        return Vrpu_tqs;
    }
}

void flexmat6::update_as_rpu_tsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNr - vu*iNr*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vq*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rpu_tsq(){
    if(Nrpu_tsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vp*iNr + vu*iNp*iNr;
        locations.col(1) = vt + vs*iNt + vq*iNs*iNt;
        Vrpu_tsq = sp_mat(locations.t(), vValues, iNr*iNp*iNu,iNt*iNs*iNq);
        Nrpu_tsq = 1;
        return Vrpu_tsq;
    }
    else{
        return Vrpu_tsq;
    }
}

void flexmat6::update_as_rqp_stu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vp*iNr*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vu*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqp_stu(){
    if(Nrqp_stu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vp*iNq*iNr;
        locations.col(1) = vs + vt*iNs + vu*iNt*iNs;
        Vrqp_stu = sp_mat(locations.t(), vValues, iNr*iNq*iNp,iNs*iNt*iNu);
        Nrqp_stu = 1;
        return Vrqp_stu;
    }
    else{
        return Vrqp_stu;
    }
}

void flexmat6::update_as_rqp_sut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vp*iNr*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vt*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqp_sut(){
    if(Nrqp_sut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vp*iNq*iNr;
        locations.col(1) = vs + vu*iNs + vt*iNu*iNs;
        Vrqp_sut = sp_mat(locations.t(), vValues, iNr*iNq*iNp,iNs*iNu*iNt);
        Nrqp_sut = 1;
        return Vrqp_sut;
    }
    else{
        return Vrqp_sut;
    }
}

void flexmat6::update_as_rqp_tsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vp*iNr*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vu*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqp_tsu(){
    if(Nrqp_tsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vp*iNq*iNr;
        locations.col(1) = vt + vs*iNt + vu*iNs*iNt;
        Vrqp_tsu = sp_mat(locations.t(), vValues, iNr*iNq*iNp,iNt*iNs*iNu);
        Nrqp_tsu = 1;
        return Vrqp_tsu;
    }
    else{
        return Vrqp_tsu;
    }
}

void flexmat6::update_as_rqp_tus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vp*iNr*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vs*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqp_tus(){
    if(Nrqp_tus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vp*iNq*iNr;
        locations.col(1) = vt + vu*iNt + vs*iNu*iNt;
        Vrqp_tus = sp_mat(locations.t(), vValues, iNr*iNq*iNp,iNt*iNu*iNs);
        Nrqp_tus = 1;
        return Vrqp_tus;
    }
    else{
        return Vrqp_tus;
    }
}

void flexmat6::update_as_rqp_ust(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vp*iNr*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vt*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqp_ust(){
    if(Nrqp_ust == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vp*iNq*iNr;
        locations.col(1) = vu + vs*iNu + vt*iNs*iNu;
        Vrqp_ust = sp_mat(locations.t(), vValues, iNr*iNq*iNp,iNu*iNs*iNt);
        Nrqp_ust = 1;
        return Vrqp_ust;
    }
    else{
        return Vrqp_ust;
    }
}

void flexmat6::update_as_rqp_uts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vp*iNr*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vs*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqp_uts(){
    if(Nrqp_uts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vp*iNq*iNr;
        locations.col(1) = vu + vt*iNu + vs*iNt*iNu;
        Vrqp_uts = sp_mat(locations.t(), vValues, iNr*iNq*iNp,iNu*iNt*iNs);
        Nrqp_uts = 1;
        return Vrqp_uts;
    }
    else{
        return Vrqp_uts;
    }
}

void flexmat6::update_as_rqs_ptu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vs*iNr*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vu*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqs_ptu(){
    if(Nrqs_ptu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vs*iNq*iNr;
        locations.col(1) = vp + vt*iNp + vu*iNt*iNp;
        Vrqs_ptu = sp_mat(locations.t(), vValues, iNr*iNq*iNs,iNp*iNt*iNu);
        Nrqs_ptu = 1;
        return Vrqs_ptu;
    }
    else{
        return Vrqs_ptu;
    }
}

void flexmat6::update_as_rqs_put(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vs*iNr*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vt*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqs_put(){
    if(Nrqs_put == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vs*iNq*iNr;
        locations.col(1) = vp + vu*iNp + vt*iNu*iNp;
        Vrqs_put = sp_mat(locations.t(), vValues, iNr*iNq*iNs,iNp*iNu*iNt);
        Nrqs_put = 1;
        return Vrqs_put;
    }
    else{
        return Vrqs_put;
    }
}

void flexmat6::update_as_rqs_tpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vs*iNr*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vu*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqs_tpu(){
    if(Nrqs_tpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vs*iNq*iNr;
        locations.col(1) = vt + vp*iNt + vu*iNp*iNt;
        Vrqs_tpu = sp_mat(locations.t(), vValues, iNr*iNq*iNs,iNt*iNp*iNu);
        Nrqs_tpu = 1;
        return Vrqs_tpu;
    }
    else{
        return Vrqs_tpu;
    }
}

void flexmat6::update_as_rqs_tup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vs*iNr*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vp*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqs_tup(){
    if(Nrqs_tup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vs*iNq*iNr;
        locations.col(1) = vt + vu*iNt + vp*iNu*iNt;
        Vrqs_tup = sp_mat(locations.t(), vValues, iNr*iNq*iNs,iNt*iNu*iNp);
        Nrqs_tup = 1;
        return Vrqs_tup;
    }
    else{
        return Vrqs_tup;
    }
}

void flexmat6::update_as_rqs_upt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vs*iNr*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vt*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqs_upt(){
    if(Nrqs_upt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vs*iNq*iNr;
        locations.col(1) = vu + vp*iNu + vt*iNp*iNu;
        Vrqs_upt = sp_mat(locations.t(), vValues, iNr*iNq*iNs,iNu*iNp*iNt);
        Nrqs_upt = 1;
        return Vrqs_upt;
    }
    else{
        return Vrqs_upt;
    }
}

void flexmat6::update_as_rqs_utp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vs*iNr*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vp*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqs_utp(){
    if(Nrqs_utp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vs*iNq*iNr;
        locations.col(1) = vu + vt*iNu + vp*iNt*iNu;
        Vrqs_utp = sp_mat(locations.t(), vValues, iNr*iNq*iNs,iNu*iNt*iNp);
        Nrqs_utp = 1;
        return Vrqs_utp;
    }
    else{
        return Vrqs_utp;
    }
}

void flexmat6::update_as_rqt_psu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vt*iNr*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vu*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqt_psu(){
    if(Nrqt_psu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vt*iNq*iNr;
        locations.col(1) = vp + vs*iNp + vu*iNs*iNp;
        Vrqt_psu = sp_mat(locations.t(), vValues, iNr*iNq*iNt,iNp*iNs*iNu);
        Nrqt_psu = 1;
        return Vrqt_psu;
    }
    else{
        return Vrqt_psu;
    }
}

void flexmat6::update_as_rqt_pus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vt*iNr*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vs*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqt_pus(){
    if(Nrqt_pus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vt*iNq*iNr;
        locations.col(1) = vp + vu*iNp + vs*iNu*iNp;
        Vrqt_pus = sp_mat(locations.t(), vValues, iNr*iNq*iNt,iNp*iNu*iNs);
        Nrqt_pus = 1;
        return Vrqt_pus;
    }
    else{
        return Vrqt_pus;
    }
}

void flexmat6::update_as_rqt_spu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vt*iNr*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vu*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqt_spu(){
    if(Nrqt_spu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vt*iNq*iNr;
        locations.col(1) = vs + vp*iNs + vu*iNp*iNs;
        Vrqt_spu = sp_mat(locations.t(), vValues, iNr*iNq*iNt,iNs*iNp*iNu);
        Nrqt_spu = 1;
        return Vrqt_spu;
    }
    else{
        return Vrqt_spu;
    }
}

void flexmat6::update_as_rqt_sup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vt*iNr*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vp*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqt_sup(){
    if(Nrqt_sup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vt*iNq*iNr;
        locations.col(1) = vs + vu*iNs + vp*iNu*iNs;
        Vrqt_sup = sp_mat(locations.t(), vValues, iNr*iNq*iNt,iNs*iNu*iNp);
        Nrqt_sup = 1;
        return Vrqt_sup;
    }
    else{
        return Vrqt_sup;
    }
}

void flexmat6::update_as_rqt_ups(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vt*iNr*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vs*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqt_ups(){
    if(Nrqt_ups == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vt*iNq*iNr;
        locations.col(1) = vu + vp*iNu + vs*iNp*iNu;
        Vrqt_ups = sp_mat(locations.t(), vValues, iNr*iNq*iNt,iNu*iNp*iNs);
        Nrqt_ups = 1;
        return Vrqt_ups;
    }
    else{
        return Vrqt_ups;
    }
}

void flexmat6::update_as_rqt_usp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vt*iNr*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vp*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqt_usp(){
    if(Nrqt_usp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vt*iNq*iNr;
        locations.col(1) = vu + vs*iNu + vp*iNs*iNu;
        Vrqt_usp = sp_mat(locations.t(), vValues, iNr*iNq*iNt,iNu*iNs*iNp);
        Nrqt_usp = 1;
        return Vrqt_usp;
    }
    else{
        return Vrqt_usp;
    }
}

void flexmat6::update_as_rqu_pst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vu*iNr*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vt*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqu_pst(){
    if(Nrqu_pst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vu*iNq*iNr;
        locations.col(1) = vp + vs*iNp + vt*iNs*iNp;
        Vrqu_pst = sp_mat(locations.t(), vValues, iNr*iNq*iNu,iNp*iNs*iNt);
        Nrqu_pst = 1;
        return Vrqu_pst;
    }
    else{
        return Vrqu_pst;
    }
}

void flexmat6::update_as_rqu_pts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vu*iNr*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vs*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqu_pts(){
    if(Nrqu_pts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vu*iNq*iNr;
        locations.col(1) = vp + vt*iNp + vs*iNt*iNp;
        Vrqu_pts = sp_mat(locations.t(), vValues, iNr*iNq*iNu,iNp*iNt*iNs);
        Nrqu_pts = 1;
        return Vrqu_pts;
    }
    else{
        return Vrqu_pts;
    }
}

void flexmat6::update_as_rqu_spt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vu*iNr*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vt*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqu_spt(){
    if(Nrqu_spt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vu*iNq*iNr;
        locations.col(1) = vs + vp*iNs + vt*iNp*iNs;
        Vrqu_spt = sp_mat(locations.t(), vValues, iNr*iNq*iNu,iNs*iNp*iNt);
        Nrqu_spt = 1;
        return Vrqu_spt;
    }
    else{
        return Vrqu_spt;
    }
}

void flexmat6::update_as_rqu_stp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vu*iNr*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vp*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqu_stp(){
    if(Nrqu_stp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vu*iNq*iNr;
        locations.col(1) = vs + vt*iNs + vp*iNt*iNs;
        Vrqu_stp = sp_mat(locations.t(), vValues, iNr*iNq*iNu,iNs*iNt*iNp);
        Nrqu_stp = 1;
        return Vrqu_stp;
    }
    else{
        return Vrqu_stp;
    }
}

void flexmat6::update_as_rqu_tps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vu*iNr*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vs*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqu_tps(){
    if(Nrqu_tps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vu*iNq*iNr;
        locations.col(1) = vt + vp*iNt + vs*iNp*iNt;
        Vrqu_tps = sp_mat(locations.t(), vValues, iNr*iNq*iNu,iNt*iNp*iNs);
        Nrqu_tps = 1;
        return Vrqu_tps;
    }
    else{
        return Vrqu_tps;
    }
}

void flexmat6::update_as_rqu_tsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNr - vu*iNr*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vp*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rqu_tsp(){
    if(Nrqu_tsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vq*iNr + vu*iNq*iNr;
        locations.col(1) = vt + vs*iNt + vp*iNs*iNt;
        Vrqu_tsp = sp_mat(locations.t(), vValues, iNr*iNq*iNu,iNt*iNs*iNp);
        Nrqu_tsp = 1;
        return Vrqu_tsp;
    }
    else{
        return Vrqu_tsp;
    }
}

void flexmat6::update_as_rsp_qtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vp*iNr*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vu*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsp_qtu(){
    if(Nrsp_qtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vp*iNs*iNr;
        locations.col(1) = vq + vt*iNq + vu*iNt*iNq;
        Vrsp_qtu = sp_mat(locations.t(), vValues, iNr*iNs*iNp,iNq*iNt*iNu);
        Nrsp_qtu = 1;
        return Vrsp_qtu;
    }
    else{
        return Vrsp_qtu;
    }
}

void flexmat6::update_as_rsp_qut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vp*iNr*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vt*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsp_qut(){
    if(Nrsp_qut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vp*iNs*iNr;
        locations.col(1) = vq + vu*iNq + vt*iNu*iNq;
        Vrsp_qut = sp_mat(locations.t(), vValues, iNr*iNs*iNp,iNq*iNu*iNt);
        Nrsp_qut = 1;
        return Vrsp_qut;
    }
    else{
        return Vrsp_qut;
    }
}

void flexmat6::update_as_rsp_tqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vp*iNr*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vu*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsp_tqu(){
    if(Nrsp_tqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vp*iNs*iNr;
        locations.col(1) = vt + vq*iNt + vu*iNq*iNt;
        Vrsp_tqu = sp_mat(locations.t(), vValues, iNr*iNs*iNp,iNt*iNq*iNu);
        Nrsp_tqu = 1;
        return Vrsp_tqu;
    }
    else{
        return Vrsp_tqu;
    }
}

void flexmat6::update_as_rsp_tuq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vp*iNr*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vq*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsp_tuq(){
    if(Nrsp_tuq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vp*iNs*iNr;
        locations.col(1) = vt + vu*iNt + vq*iNu*iNt;
        Vrsp_tuq = sp_mat(locations.t(), vValues, iNr*iNs*iNp,iNt*iNu*iNq);
        Nrsp_tuq = 1;
        return Vrsp_tuq;
    }
    else{
        return Vrsp_tuq;
    }
}

void flexmat6::update_as_rsp_uqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vp*iNr*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vt*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsp_uqt(){
    if(Nrsp_uqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vp*iNs*iNr;
        locations.col(1) = vu + vq*iNu + vt*iNq*iNu;
        Vrsp_uqt = sp_mat(locations.t(), vValues, iNr*iNs*iNp,iNu*iNq*iNt);
        Nrsp_uqt = 1;
        return Vrsp_uqt;
    }
    else{
        return Vrsp_uqt;
    }
}

void flexmat6::update_as_rsp_utq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vp*iNr*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vq*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsp_utq(){
    if(Nrsp_utq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vp*iNs*iNr;
        locations.col(1) = vu + vt*iNu + vq*iNt*iNu;
        Vrsp_utq = sp_mat(locations.t(), vValues, iNr*iNs*iNp,iNu*iNt*iNq);
        Nrsp_utq = 1;
        return Vrsp_utq;
    }
    else{
        return Vrsp_utq;
    }
}

void flexmat6::update_as_rsq_ptu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vq*iNr*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vu*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsq_ptu(){
    if(Nrsq_ptu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vq*iNs*iNr;
        locations.col(1) = vp + vt*iNp + vu*iNt*iNp;
        Vrsq_ptu = sp_mat(locations.t(), vValues, iNr*iNs*iNq,iNp*iNt*iNu);
        Nrsq_ptu = 1;
        return Vrsq_ptu;
    }
    else{
        return Vrsq_ptu;
    }
}

void flexmat6::update_as_rsq_put(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vq*iNr*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vt*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsq_put(){
    if(Nrsq_put == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vq*iNs*iNr;
        locations.col(1) = vp + vu*iNp + vt*iNu*iNp;
        Vrsq_put = sp_mat(locations.t(), vValues, iNr*iNs*iNq,iNp*iNu*iNt);
        Nrsq_put = 1;
        return Vrsq_put;
    }
    else{
        return Vrsq_put;
    }
}

void flexmat6::update_as_rsq_tpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vq*iNr*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vu*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsq_tpu(){
    if(Nrsq_tpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vq*iNs*iNr;
        locations.col(1) = vt + vp*iNt + vu*iNp*iNt;
        Vrsq_tpu = sp_mat(locations.t(), vValues, iNr*iNs*iNq,iNt*iNp*iNu);
        Nrsq_tpu = 1;
        return Vrsq_tpu;
    }
    else{
        return Vrsq_tpu;
    }
}

void flexmat6::update_as_rsq_tup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vq*iNr*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vp*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsq_tup(){
    if(Nrsq_tup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vq*iNs*iNr;
        locations.col(1) = vt + vu*iNt + vp*iNu*iNt;
        Vrsq_tup = sp_mat(locations.t(), vValues, iNr*iNs*iNq,iNt*iNu*iNp);
        Nrsq_tup = 1;
        return Vrsq_tup;
    }
    else{
        return Vrsq_tup;
    }
}

void flexmat6::update_as_rsq_upt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vq*iNr*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vt*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsq_upt(){
    if(Nrsq_upt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vq*iNs*iNr;
        locations.col(1) = vu + vp*iNu + vt*iNp*iNu;
        Vrsq_upt = sp_mat(locations.t(), vValues, iNr*iNs*iNq,iNu*iNp*iNt);
        Nrsq_upt = 1;
        return Vrsq_upt;
    }
    else{
        return Vrsq_upt;
    }
}

void flexmat6::update_as_rsq_utp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vq*iNr*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vp*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsq_utp(){
    if(Nrsq_utp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vq*iNs*iNr;
        locations.col(1) = vu + vt*iNu + vp*iNt*iNu;
        Vrsq_utp = sp_mat(locations.t(), vValues, iNr*iNs*iNq,iNu*iNt*iNp);
        Nrsq_utp = 1;
        return Vrsq_utp;
    }
    else{
        return Vrsq_utp;
    }
}

void flexmat6::update_as_rst_pqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vt*iNr*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vu*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rst_pqu(){
    if(Nrst_pqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vt*iNs*iNr;
        locations.col(1) = vp + vq*iNp + vu*iNq*iNp;
        Vrst_pqu = sp_mat(locations.t(), vValues, iNr*iNs*iNt,iNp*iNq*iNu);
        Nrst_pqu = 1;
        return Vrst_pqu;
    }
    else{
        return Vrst_pqu;
    }
}

void flexmat6::update_as_rst_puq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vt*iNr*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vq*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rst_puq(){
    if(Nrst_puq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vt*iNs*iNr;
        locations.col(1) = vp + vu*iNp + vq*iNu*iNp;
        Vrst_puq = sp_mat(locations.t(), vValues, iNr*iNs*iNt,iNp*iNu*iNq);
        Nrst_puq = 1;
        return Vrst_puq;
    }
    else{
        return Vrst_puq;
    }
}

void flexmat6::update_as_rst_qpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vt*iNr*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vu*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rst_qpu(){
    if(Nrst_qpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vt*iNs*iNr;
        locations.col(1) = vq + vp*iNq + vu*iNp*iNq;
        Vrst_qpu = sp_mat(locations.t(), vValues, iNr*iNs*iNt,iNq*iNp*iNu);
        Nrst_qpu = 1;
        return Vrst_qpu;
    }
    else{
        return Vrst_qpu;
    }
}

void flexmat6::update_as_rst_qup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vt*iNr*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vp*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rst_qup(){
    if(Nrst_qup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vt*iNs*iNr;
        locations.col(1) = vq + vu*iNq + vp*iNu*iNq;
        Vrst_qup = sp_mat(locations.t(), vValues, iNr*iNs*iNt,iNq*iNu*iNp);
        Nrst_qup = 1;
        return Vrst_qup;
    }
    else{
        return Vrst_qup;
    }
}

void flexmat6::update_as_rst_upq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vt*iNr*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vq*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rst_upq(){
    if(Nrst_upq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vt*iNs*iNr;
        locations.col(1) = vu + vp*iNu + vq*iNp*iNu;
        Vrst_upq = sp_mat(locations.t(), vValues, iNr*iNs*iNt,iNu*iNp*iNq);
        Nrst_upq = 1;
        return Vrst_upq;
    }
    else{
        return Vrst_upq;
    }
}

void flexmat6::update_as_rst_uqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vt*iNr*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vp*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rst_uqp(){
    if(Nrst_uqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vt*iNs*iNr;
        locations.col(1) = vu + vq*iNu + vp*iNq*iNu;
        Vrst_uqp = sp_mat(locations.t(), vValues, iNr*iNs*iNt,iNu*iNq*iNp);
        Nrst_uqp = 1;
        return Vrst_uqp;
    }
    else{
        return Vrst_uqp;
    }
}

void flexmat6::update_as_rsu_pqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vu*iNr*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vt*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsu_pqt(){
    if(Nrsu_pqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vu*iNs*iNr;
        locations.col(1) = vp + vq*iNp + vt*iNq*iNp;
        Vrsu_pqt = sp_mat(locations.t(), vValues, iNr*iNs*iNu,iNp*iNq*iNt);
        Nrsu_pqt = 1;
        return Vrsu_pqt;
    }
    else{
        return Vrsu_pqt;
    }
}

void flexmat6::update_as_rsu_ptq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vu*iNr*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vq*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsu_ptq(){
    if(Nrsu_ptq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vu*iNs*iNr;
        locations.col(1) = vp + vt*iNp + vq*iNt*iNp;
        Vrsu_ptq = sp_mat(locations.t(), vValues, iNr*iNs*iNu,iNp*iNt*iNq);
        Nrsu_ptq = 1;
        return Vrsu_ptq;
    }
    else{
        return Vrsu_ptq;
    }
}

void flexmat6::update_as_rsu_qpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vu*iNr*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vt*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsu_qpt(){
    if(Nrsu_qpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vu*iNs*iNr;
        locations.col(1) = vq + vp*iNq + vt*iNp*iNq;
        Vrsu_qpt = sp_mat(locations.t(), vValues, iNr*iNs*iNu,iNq*iNp*iNt);
        Nrsu_qpt = 1;
        return Vrsu_qpt;
    }
    else{
        return Vrsu_qpt;
    }
}

void flexmat6::update_as_rsu_qtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vu*iNr*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vp*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsu_qtp(){
    if(Nrsu_qtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vu*iNs*iNr;
        locations.col(1) = vq + vt*iNq + vp*iNt*iNq;
        Vrsu_qtp = sp_mat(locations.t(), vValues, iNr*iNs*iNu,iNq*iNt*iNp);
        Nrsu_qtp = 1;
        return Vrsu_qtp;
    }
    else{
        return Vrsu_qtp;
    }
}

void flexmat6::update_as_rsu_tpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vu*iNr*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vq*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsu_tpq(){
    if(Nrsu_tpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vu*iNs*iNr;
        locations.col(1) = vt + vp*iNt + vq*iNp*iNt;
        Vrsu_tpq = sp_mat(locations.t(), vValues, iNr*iNs*iNu,iNt*iNp*iNq);
        Nrsu_tpq = 1;
        return Vrsu_tpq;
    }
    else{
        return Vrsu_tpq;
    }
}

void flexmat6::update_as_rsu_tqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNr - vu*iNr*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vp*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rsu_tqp(){
    if(Nrsu_tqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vs*iNr + vu*iNs*iNr;
        locations.col(1) = vt + vq*iNt + vp*iNq*iNt;
        Vrsu_tqp = sp_mat(locations.t(), vValues, iNr*iNs*iNu,iNt*iNq*iNp);
        Nrsu_tqp = 1;
        return Vrsu_tqp;
    }
    else{
        return Vrsu_tqp;
    }
}

void flexmat6::update_as_rtp_qsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vp*iNr*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vu*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtp_qsu(){
    if(Nrtp_qsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vp*iNt*iNr;
        locations.col(1) = vq + vs*iNq + vu*iNs*iNq;
        Vrtp_qsu = sp_mat(locations.t(), vValues, iNr*iNt*iNp,iNq*iNs*iNu);
        Nrtp_qsu = 1;
        return Vrtp_qsu;
    }
    else{
        return Vrtp_qsu;
    }
}

void flexmat6::update_as_rtp_qus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vp*iNr*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vs*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtp_qus(){
    if(Nrtp_qus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vp*iNt*iNr;
        locations.col(1) = vq + vu*iNq + vs*iNu*iNq;
        Vrtp_qus = sp_mat(locations.t(), vValues, iNr*iNt*iNp,iNq*iNu*iNs);
        Nrtp_qus = 1;
        return Vrtp_qus;
    }
    else{
        return Vrtp_qus;
    }
}

void flexmat6::update_as_rtp_squ(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vp*iNr*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vu*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtp_squ(){
    if(Nrtp_squ == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vp*iNt*iNr;
        locations.col(1) = vs + vq*iNs + vu*iNq*iNs;
        Vrtp_squ = sp_mat(locations.t(), vValues, iNr*iNt*iNp,iNs*iNq*iNu);
        Nrtp_squ = 1;
        return Vrtp_squ;
    }
    else{
        return Vrtp_squ;
    }
}

void flexmat6::update_as_rtp_suq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vp*iNr*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vq*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtp_suq(){
    if(Nrtp_suq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vp*iNt*iNr;
        locations.col(1) = vs + vu*iNs + vq*iNu*iNs;
        Vrtp_suq = sp_mat(locations.t(), vValues, iNr*iNt*iNp,iNs*iNu*iNq);
        Nrtp_suq = 1;
        return Vrtp_suq;
    }
    else{
        return Vrtp_suq;
    }
}

void flexmat6::update_as_rtp_uqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vp*iNr*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vs*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtp_uqs(){
    if(Nrtp_uqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vp*iNt*iNr;
        locations.col(1) = vu + vq*iNu + vs*iNq*iNu;
        Vrtp_uqs = sp_mat(locations.t(), vValues, iNr*iNt*iNp,iNu*iNq*iNs);
        Nrtp_uqs = 1;
        return Vrtp_uqs;
    }
    else{
        return Vrtp_uqs;
    }
}

void flexmat6::update_as_rtp_usq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vp*iNr*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vq*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtp_usq(){
    if(Nrtp_usq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vp*iNt*iNr;
        locations.col(1) = vu + vs*iNu + vq*iNs*iNu;
        Vrtp_usq = sp_mat(locations.t(), vValues, iNr*iNt*iNp,iNu*iNs*iNq);
        Nrtp_usq = 1;
        return Vrtp_usq;
    }
    else{
        return Vrtp_usq;
    }
}

void flexmat6::update_as_rtq_psu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vq*iNr*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vu*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtq_psu(){
    if(Nrtq_psu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vq*iNt*iNr;
        locations.col(1) = vp + vs*iNp + vu*iNs*iNp;
        Vrtq_psu = sp_mat(locations.t(), vValues, iNr*iNt*iNq,iNp*iNs*iNu);
        Nrtq_psu = 1;
        return Vrtq_psu;
    }
    else{
        return Vrtq_psu;
    }
}

void flexmat6::update_as_rtq_pus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vq*iNr*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vs*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtq_pus(){
    if(Nrtq_pus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vq*iNt*iNr;
        locations.col(1) = vp + vu*iNp + vs*iNu*iNp;
        Vrtq_pus = sp_mat(locations.t(), vValues, iNr*iNt*iNq,iNp*iNu*iNs);
        Nrtq_pus = 1;
        return Vrtq_pus;
    }
    else{
        return Vrtq_pus;
    }
}

void flexmat6::update_as_rtq_spu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vq*iNr*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vu*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtq_spu(){
    if(Nrtq_spu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vq*iNt*iNr;
        locations.col(1) = vs + vp*iNs + vu*iNp*iNs;
        Vrtq_spu = sp_mat(locations.t(), vValues, iNr*iNt*iNq,iNs*iNp*iNu);
        Nrtq_spu = 1;
        return Vrtq_spu;
    }
    else{
        return Vrtq_spu;
    }
}

void flexmat6::update_as_rtq_sup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vq*iNr*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vp*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtq_sup(){
    if(Nrtq_sup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vq*iNt*iNr;
        locations.col(1) = vs + vu*iNs + vp*iNu*iNs;
        Vrtq_sup = sp_mat(locations.t(), vValues, iNr*iNt*iNq,iNs*iNu*iNp);
        Nrtq_sup = 1;
        return Vrtq_sup;
    }
    else{
        return Vrtq_sup;
    }
}

void flexmat6::update_as_rtq_ups(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vq*iNr*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vs*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtq_ups(){
    if(Nrtq_ups == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vq*iNt*iNr;
        locations.col(1) = vu + vp*iNu + vs*iNp*iNu;
        Vrtq_ups = sp_mat(locations.t(), vValues, iNr*iNt*iNq,iNu*iNp*iNs);
        Nrtq_ups = 1;
        return Vrtq_ups;
    }
    else{
        return Vrtq_ups;
    }
}

void flexmat6::update_as_rtq_usp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vq*iNr*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vp*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtq_usp(){
    if(Nrtq_usp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vq*iNt*iNr;
        locations.col(1) = vu + vs*iNu + vp*iNs*iNu;
        Vrtq_usp = sp_mat(locations.t(), vValues, iNr*iNt*iNq,iNu*iNs*iNp);
        Nrtq_usp = 1;
        return Vrtq_usp;
    }
    else{
        return Vrtq_usp;
    }
}

void flexmat6::update_as_rts_pqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vs*iNr*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vu*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rts_pqu(){
    if(Nrts_pqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vs*iNt*iNr;
        locations.col(1) = vp + vq*iNp + vu*iNq*iNp;
        Vrts_pqu = sp_mat(locations.t(), vValues, iNr*iNt*iNs,iNp*iNq*iNu);
        Nrts_pqu = 1;
        return Vrts_pqu;
    }
    else{
        return Vrts_pqu;
    }
}

void flexmat6::update_as_rts_puq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vs*iNr*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vq*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rts_puq(){
    if(Nrts_puq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vs*iNt*iNr;
        locations.col(1) = vp + vu*iNp + vq*iNu*iNp;
        Vrts_puq = sp_mat(locations.t(), vValues, iNr*iNt*iNs,iNp*iNu*iNq);
        Nrts_puq = 1;
        return Vrts_puq;
    }
    else{
        return Vrts_puq;
    }
}

void flexmat6::update_as_rts_qpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vs*iNr*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vu*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rts_qpu(){
    if(Nrts_qpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vs*iNt*iNr;
        locations.col(1) = vq + vp*iNq + vu*iNp*iNq;
        Vrts_qpu = sp_mat(locations.t(), vValues, iNr*iNt*iNs,iNq*iNp*iNu);
        Nrts_qpu = 1;
        return Vrts_qpu;
    }
    else{
        return Vrts_qpu;
    }
}

void flexmat6::update_as_rts_qup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vs*iNr*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vp*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rts_qup(){
    if(Nrts_qup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vs*iNt*iNr;
        locations.col(1) = vq + vu*iNq + vp*iNu*iNq;
        Vrts_qup = sp_mat(locations.t(), vValues, iNr*iNt*iNs,iNq*iNu*iNp);
        Nrts_qup = 1;
        return Vrts_qup;
    }
    else{
        return Vrts_qup;
    }
}

void flexmat6::update_as_rts_upq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vs*iNr*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vq*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rts_upq(){
    if(Nrts_upq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vs*iNt*iNr;
        locations.col(1) = vu + vp*iNu + vq*iNp*iNu;
        Vrts_upq = sp_mat(locations.t(), vValues, iNr*iNt*iNs,iNu*iNp*iNq);
        Nrts_upq = 1;
        return Vrts_upq;
    }
    else{
        return Vrts_upq;
    }
}

void flexmat6::update_as_rts_uqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vs*iNr*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vp*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rts_uqp(){
    if(Nrts_uqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vs*iNt*iNr;
        locations.col(1) = vu + vq*iNu + vp*iNq*iNu;
        Vrts_uqp = sp_mat(locations.t(), vValues, iNr*iNt*iNs,iNu*iNq*iNp);
        Nrts_uqp = 1;
        return Vrts_uqp;
    }
    else{
        return Vrts_uqp;
    }
}

void flexmat6::update_as_rtu_pqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vu*iNr*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vs*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtu_pqs(){
    if(Nrtu_pqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vu*iNt*iNr;
        locations.col(1) = vp + vq*iNp + vs*iNq*iNp;
        Vrtu_pqs = sp_mat(locations.t(), vValues, iNr*iNt*iNu,iNp*iNq*iNs);
        Nrtu_pqs = 1;
        return Vrtu_pqs;
    }
    else{
        return Vrtu_pqs;
    }
}

void flexmat6::update_as_rtu_psq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vu*iNr*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vq*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtu_psq(){
    if(Nrtu_psq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vu*iNt*iNr;
        locations.col(1) = vp + vs*iNp + vq*iNs*iNp;
        Vrtu_psq = sp_mat(locations.t(), vValues, iNr*iNt*iNu,iNp*iNs*iNq);
        Nrtu_psq = 1;
        return Vrtu_psq;
    }
    else{
        return Vrtu_psq;
    }
}

void flexmat6::update_as_rtu_qps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vu*iNr*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vs*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtu_qps(){
    if(Nrtu_qps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vu*iNt*iNr;
        locations.col(1) = vq + vp*iNq + vs*iNp*iNq;
        Vrtu_qps = sp_mat(locations.t(), vValues, iNr*iNt*iNu,iNq*iNp*iNs);
        Nrtu_qps = 1;
        return Vrtu_qps;
    }
    else{
        return Vrtu_qps;
    }
}

void flexmat6::update_as_rtu_qsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vu*iNr*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vp*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtu_qsp(){
    if(Nrtu_qsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vu*iNt*iNr;
        locations.col(1) = vq + vs*iNq + vp*iNs*iNq;
        Vrtu_qsp = sp_mat(locations.t(), vValues, iNr*iNt*iNu,iNq*iNs*iNp);
        Nrtu_qsp = 1;
        return Vrtu_qsp;
    }
    else{
        return Vrtu_qsp;
    }
}

void flexmat6::update_as_rtu_spq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vu*iNr*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vq*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtu_spq(){
    if(Nrtu_spq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vu*iNt*iNr;
        locations.col(1) = vs + vp*iNs + vq*iNp*iNs;
        Vrtu_spq = sp_mat(locations.t(), vValues, iNr*iNt*iNu,iNs*iNp*iNq);
        Nrtu_spq = 1;
        return Vrtu_spq;
    }
    else{
        return Vrtu_spq;
    }
}

void flexmat6::update_as_rtu_sqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNr - vu*iNr*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vp*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rtu_sqp(){
    if(Nrtu_sqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vt*iNr + vu*iNt*iNr;
        locations.col(1) = vs + vq*iNs + vp*iNq*iNs;
        Vrtu_sqp = sp_mat(locations.t(), vValues, iNr*iNt*iNu,iNs*iNq*iNp);
        Nrtu_sqp = 1;
        return Vrtu_sqp;
    }
    else{
        return Vrtu_sqp;
    }
}

void flexmat6::update_as_rup_qst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vp*iNr*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vt*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rup_qst(){
    if(Nrup_qst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vp*iNu*iNr;
        locations.col(1) = vq + vs*iNq + vt*iNs*iNq;
        Vrup_qst = sp_mat(locations.t(), vValues, iNr*iNu*iNp,iNq*iNs*iNt);
        Nrup_qst = 1;
        return Vrup_qst;
    }
    else{
        return Vrup_qst;
    }
}

void flexmat6::update_as_rup_qts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vp*iNr*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vs*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rup_qts(){
    if(Nrup_qts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vp*iNu*iNr;
        locations.col(1) = vq + vt*iNq + vs*iNt*iNq;
        Vrup_qts = sp_mat(locations.t(), vValues, iNr*iNu*iNp,iNq*iNt*iNs);
        Nrup_qts = 1;
        return Vrup_qts;
    }
    else{
        return Vrup_qts;
    }
}

void flexmat6::update_as_rup_sqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vp*iNr*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vt*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rup_sqt(){
    if(Nrup_sqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vp*iNu*iNr;
        locations.col(1) = vs + vq*iNs + vt*iNq*iNs;
        Vrup_sqt = sp_mat(locations.t(), vValues, iNr*iNu*iNp,iNs*iNq*iNt);
        Nrup_sqt = 1;
        return Vrup_sqt;
    }
    else{
        return Vrup_sqt;
    }
}

void flexmat6::update_as_rup_stq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vp*iNr*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vq*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rup_stq(){
    if(Nrup_stq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vp*iNu*iNr;
        locations.col(1) = vs + vt*iNs + vq*iNt*iNs;
        Vrup_stq = sp_mat(locations.t(), vValues, iNr*iNu*iNp,iNs*iNt*iNq);
        Nrup_stq = 1;
        return Vrup_stq;
    }
    else{
        return Vrup_stq;
    }
}

void flexmat6::update_as_rup_tqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vp*iNr*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vs*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rup_tqs(){
    if(Nrup_tqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vp*iNu*iNr;
        locations.col(1) = vt + vq*iNt + vs*iNq*iNt;
        Vrup_tqs = sp_mat(locations.t(), vValues, iNr*iNu*iNp,iNt*iNq*iNs);
        Nrup_tqs = 1;
        return Vrup_tqs;
    }
    else{
        return Vrup_tqs;
    }
}

void flexmat6::update_as_rup_tsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vp*iNr*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vq*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rup_tsq(){
    if(Nrup_tsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vp*iNu*iNr;
        locations.col(1) = vt + vs*iNt + vq*iNs*iNt;
        Vrup_tsq = sp_mat(locations.t(), vValues, iNr*iNu*iNp,iNt*iNs*iNq);
        Nrup_tsq = 1;
        return Vrup_tsq;
    }
    else{
        return Vrup_tsq;
    }
}

void flexmat6::update_as_ruq_pst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vq*iNr*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vt*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ruq_pst(){
    if(Nruq_pst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vq*iNu*iNr;
        locations.col(1) = vp + vs*iNp + vt*iNs*iNp;
        Vruq_pst = sp_mat(locations.t(), vValues, iNr*iNu*iNq,iNp*iNs*iNt);
        Nruq_pst = 1;
        return Vruq_pst;
    }
    else{
        return Vruq_pst;
    }
}

void flexmat6::update_as_ruq_pts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vq*iNr*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vs*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ruq_pts(){
    if(Nruq_pts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vq*iNu*iNr;
        locations.col(1) = vp + vt*iNp + vs*iNt*iNp;
        Vruq_pts = sp_mat(locations.t(), vValues, iNr*iNu*iNq,iNp*iNt*iNs);
        Nruq_pts = 1;
        return Vruq_pts;
    }
    else{
        return Vruq_pts;
    }
}

void flexmat6::update_as_ruq_spt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vq*iNr*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vt*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ruq_spt(){
    if(Nruq_spt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vq*iNu*iNr;
        locations.col(1) = vs + vp*iNs + vt*iNp*iNs;
        Vruq_spt = sp_mat(locations.t(), vValues, iNr*iNu*iNq,iNs*iNp*iNt);
        Nruq_spt = 1;
        return Vruq_spt;
    }
    else{
        return Vruq_spt;
    }
}

void flexmat6::update_as_ruq_stp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vq*iNr*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vp*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ruq_stp(){
    if(Nruq_stp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vq*iNu*iNr;
        locations.col(1) = vs + vt*iNs + vp*iNt*iNs;
        Vruq_stp = sp_mat(locations.t(), vValues, iNr*iNu*iNq,iNs*iNt*iNp);
        Nruq_stp = 1;
        return Vruq_stp;
    }
    else{
        return Vruq_stp;
    }
}

void flexmat6::update_as_ruq_tps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vq*iNr*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vs*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ruq_tps(){
    if(Nruq_tps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vq*iNu*iNr;
        locations.col(1) = vt + vp*iNt + vs*iNp*iNt;
        Vruq_tps = sp_mat(locations.t(), vValues, iNr*iNu*iNq,iNt*iNp*iNs);
        Nruq_tps = 1;
        return Vruq_tps;
    }
    else{
        return Vruq_tps;
    }
}

void flexmat6::update_as_ruq_tsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vq*iNr*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vp*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ruq_tsp(){
    if(Nruq_tsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vq*iNu*iNr;
        locations.col(1) = vt + vs*iNt + vp*iNs*iNt;
        Vruq_tsp = sp_mat(locations.t(), vValues, iNr*iNu*iNq,iNt*iNs*iNp);
        Nruq_tsp = 1;
        return Vruq_tsp;
    }
    else{
        return Vruq_tsp;
    }
}

void flexmat6::update_as_rus_pqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vs*iNr*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vt*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rus_pqt(){
    if(Nrus_pqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vs*iNu*iNr;
        locations.col(1) = vp + vq*iNp + vt*iNq*iNp;
        Vrus_pqt = sp_mat(locations.t(), vValues, iNr*iNu*iNs,iNp*iNq*iNt);
        Nrus_pqt = 1;
        return Vrus_pqt;
    }
    else{
        return Vrus_pqt;
    }
}

void flexmat6::update_as_rus_ptq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vs*iNr*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vq*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rus_ptq(){
    if(Nrus_ptq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vs*iNu*iNr;
        locations.col(1) = vp + vt*iNp + vq*iNt*iNp;
        Vrus_ptq = sp_mat(locations.t(), vValues, iNr*iNu*iNs,iNp*iNt*iNq);
        Nrus_ptq = 1;
        return Vrus_ptq;
    }
    else{
        return Vrus_ptq;
    }
}

void flexmat6::update_as_rus_qpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vs*iNr*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vt*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rus_qpt(){
    if(Nrus_qpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vs*iNu*iNr;
        locations.col(1) = vq + vp*iNq + vt*iNp*iNq;
        Vrus_qpt = sp_mat(locations.t(), vValues, iNr*iNu*iNs,iNq*iNp*iNt);
        Nrus_qpt = 1;
        return Vrus_qpt;
    }
    else{
        return Vrus_qpt;
    }
}

void flexmat6::update_as_rus_qtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vs*iNr*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vp*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rus_qtp(){
    if(Nrus_qtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vs*iNu*iNr;
        locations.col(1) = vq + vt*iNq + vp*iNt*iNq;
        Vrus_qtp = sp_mat(locations.t(), vValues, iNr*iNu*iNs,iNq*iNt*iNp);
        Nrus_qtp = 1;
        return Vrus_qtp;
    }
    else{
        return Vrus_qtp;
    }
}

void flexmat6::update_as_rus_tpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vs*iNr*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vq*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rus_tpq(){
    if(Nrus_tpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vs*iNu*iNr;
        locations.col(1) = vt + vp*iNt + vq*iNp*iNt;
        Vrus_tpq = sp_mat(locations.t(), vValues, iNr*iNu*iNs,iNt*iNp*iNq);
        Nrus_tpq = 1;
        return Vrus_tpq;
    }
    else{
        return Vrus_tpq;
    }
}

void flexmat6::update_as_rus_tqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vs*iNr*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vp*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rus_tqp(){
    if(Nrus_tqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vs*iNu*iNr;
        locations.col(1) = vt + vq*iNt + vp*iNq*iNt;
        Vrus_tqp = sp_mat(locations.t(), vValues, iNr*iNu*iNs,iNt*iNq*iNp);
        Nrus_tqp = 1;
        return Vrus_tqp;
    }
    else{
        return Vrus_tqp;
    }
}

void flexmat6::update_as_rut_pqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vt*iNr*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vs*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rut_pqs(){
    if(Nrut_pqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vt*iNu*iNr;
        locations.col(1) = vp + vq*iNp + vs*iNq*iNp;
        Vrut_pqs = sp_mat(locations.t(), vValues, iNr*iNu*iNt,iNp*iNq*iNs);
        Nrut_pqs = 1;
        return Vrut_pqs;
    }
    else{
        return Vrut_pqs;
    }
}

void flexmat6::update_as_rut_psq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vt*iNr*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vq*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rut_psq(){
    if(Nrut_psq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vt*iNu*iNr;
        locations.col(1) = vp + vs*iNp + vq*iNs*iNp;
        Vrut_psq = sp_mat(locations.t(), vValues, iNr*iNu*iNt,iNp*iNs*iNq);
        Nrut_psq = 1;
        return Vrut_psq;
    }
    else{
        return Vrut_psq;
    }
}

void flexmat6::update_as_rut_qps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vt*iNr*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vs*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rut_qps(){
    if(Nrut_qps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vt*iNu*iNr;
        locations.col(1) = vq + vp*iNq + vs*iNp*iNq;
        Vrut_qps = sp_mat(locations.t(), vValues, iNr*iNu*iNt,iNq*iNp*iNs);
        Nrut_qps = 1;
        return Vrut_qps;
    }
    else{
        return Vrut_qps;
    }
}

void flexmat6::update_as_rut_qsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vt*iNr*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vp*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rut_qsp(){
    if(Nrut_qsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vt*iNu*iNr;
        locations.col(1) = vq + vs*iNq + vp*iNs*iNq;
        Vrut_qsp = sp_mat(locations.t(), vValues, iNr*iNu*iNt,iNq*iNs*iNp);
        Nrut_qsp = 1;
        return Vrut_qsp;
    }
    else{
        return Vrut_qsp;
    }
}

void flexmat6::update_as_rut_spq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vt*iNr*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vq*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rut_spq(){
    if(Nrut_spq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vt*iNu*iNr;
        locations.col(1) = vs + vp*iNs + vq*iNp*iNs;
        Vrut_spq = sp_mat(locations.t(), vValues, iNr*iNu*iNt,iNs*iNp*iNq);
        Nrut_spq = 1;
        return Vrut_spq;
    }
    else{
        return Vrut_spq;
    }
}

void flexmat6::update_as_rut_sqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNr - vt*iNr*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vp*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::rut_sqp(){
    if(Nrut_sqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vr + vu*iNr + vt*iNu*iNr;
        locations.col(1) = vs + vq*iNs + vp*iNq*iNs;
        Vrut_sqp = sp_mat(locations.t(), vValues, iNr*iNu*iNt,iNs*iNq*iNp);
        Nrut_sqp = 1;
        return Vrut_sqp;
    }
    else{
        return Vrut_sqp;
    }
}

void flexmat6::update_as_spq_rtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vq*iNs*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vu*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spq_rtu(){
    if(Nspq_rtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vq*iNp*iNs;
        locations.col(1) = vr + vt*iNr + vu*iNt*iNr;
        Vspq_rtu = sp_mat(locations.t(), vValues, iNs*iNp*iNq,iNr*iNt*iNu);
        Nspq_rtu = 1;
        return Vspq_rtu;
    }
    else{
        return Vspq_rtu;
    }
}

void flexmat6::update_as_spq_rut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vq*iNs*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vt*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spq_rut(){
    if(Nspq_rut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vq*iNp*iNs;
        locations.col(1) = vr + vu*iNr + vt*iNu*iNr;
        Vspq_rut = sp_mat(locations.t(), vValues, iNs*iNp*iNq,iNr*iNu*iNt);
        Nspq_rut = 1;
        return Vspq_rut;
    }
    else{
        return Vspq_rut;
    }
}

void flexmat6::update_as_spq_tru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vq*iNs*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vu*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spq_tru(){
    if(Nspq_tru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vq*iNp*iNs;
        locations.col(1) = vt + vr*iNt + vu*iNr*iNt;
        Vspq_tru = sp_mat(locations.t(), vValues, iNs*iNp*iNq,iNt*iNr*iNu);
        Nspq_tru = 1;
        return Vspq_tru;
    }
    else{
        return Vspq_tru;
    }
}

void flexmat6::update_as_spq_tur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vq*iNs*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vr*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spq_tur(){
    if(Nspq_tur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vq*iNp*iNs;
        locations.col(1) = vt + vu*iNt + vr*iNu*iNt;
        Vspq_tur = sp_mat(locations.t(), vValues, iNs*iNp*iNq,iNt*iNu*iNr);
        Nspq_tur = 1;
        return Vspq_tur;
    }
    else{
        return Vspq_tur;
    }
}

void flexmat6::update_as_spq_urt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vq*iNs*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vt*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spq_urt(){
    if(Nspq_urt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vq*iNp*iNs;
        locations.col(1) = vu + vr*iNu + vt*iNr*iNu;
        Vspq_urt = sp_mat(locations.t(), vValues, iNs*iNp*iNq,iNu*iNr*iNt);
        Nspq_urt = 1;
        return Vspq_urt;
    }
    else{
        return Vspq_urt;
    }
}

void flexmat6::update_as_spq_utr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vq*iNs*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vr*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spq_utr(){
    if(Nspq_utr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vq*iNp*iNs;
        locations.col(1) = vu + vt*iNu + vr*iNt*iNu;
        Vspq_utr = sp_mat(locations.t(), vValues, iNs*iNp*iNq,iNu*iNt*iNr);
        Nspq_utr = 1;
        return Vspq_utr;
    }
    else{
        return Vspq_utr;
    }
}

void flexmat6::update_as_spr_qtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vr*iNs*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vu*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spr_qtu(){
    if(Nspr_qtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vr*iNp*iNs;
        locations.col(1) = vq + vt*iNq + vu*iNt*iNq;
        Vspr_qtu = sp_mat(locations.t(), vValues, iNs*iNp*iNr,iNq*iNt*iNu);
        Nspr_qtu = 1;
        return Vspr_qtu;
    }
    else{
        return Vspr_qtu;
    }
}

void flexmat6::update_as_spr_qut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vr*iNs*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vt*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spr_qut(){
    if(Nspr_qut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vr*iNp*iNs;
        locations.col(1) = vq + vu*iNq + vt*iNu*iNq;
        Vspr_qut = sp_mat(locations.t(), vValues, iNs*iNp*iNr,iNq*iNu*iNt);
        Nspr_qut = 1;
        return Vspr_qut;
    }
    else{
        return Vspr_qut;
    }
}

void flexmat6::update_as_spr_tqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vr*iNs*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vu*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spr_tqu(){
    if(Nspr_tqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vr*iNp*iNs;
        locations.col(1) = vt + vq*iNt + vu*iNq*iNt;
        Vspr_tqu = sp_mat(locations.t(), vValues, iNs*iNp*iNr,iNt*iNq*iNu);
        Nspr_tqu = 1;
        return Vspr_tqu;
    }
    else{
        return Vspr_tqu;
    }
}

void flexmat6::update_as_spr_tuq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vr*iNs*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vq*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spr_tuq(){
    if(Nspr_tuq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vr*iNp*iNs;
        locations.col(1) = vt + vu*iNt + vq*iNu*iNt;
        Vspr_tuq = sp_mat(locations.t(), vValues, iNs*iNp*iNr,iNt*iNu*iNq);
        Nspr_tuq = 1;
        return Vspr_tuq;
    }
    else{
        return Vspr_tuq;
    }
}

void flexmat6::update_as_spr_uqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vr*iNs*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vt*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spr_uqt(){
    if(Nspr_uqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vr*iNp*iNs;
        locations.col(1) = vu + vq*iNu + vt*iNq*iNu;
        Vspr_uqt = sp_mat(locations.t(), vValues, iNs*iNp*iNr,iNu*iNq*iNt);
        Nspr_uqt = 1;
        return Vspr_uqt;
    }
    else{
        return Vspr_uqt;
    }
}

void flexmat6::update_as_spr_utq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vr*iNs*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vq*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spr_utq(){
    if(Nspr_utq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vr*iNp*iNs;
        locations.col(1) = vu + vt*iNu + vq*iNt*iNu;
        Vspr_utq = sp_mat(locations.t(), vValues, iNs*iNp*iNr,iNu*iNt*iNq);
        Nspr_utq = 1;
        return Vspr_utq;
    }
    else{
        return Vspr_utq;
    }
}

void flexmat6::update_as_spt_qru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vt*iNs*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vu*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spt_qru(){
    if(Nspt_qru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vt*iNp*iNs;
        locations.col(1) = vq + vr*iNq + vu*iNr*iNq;
        Vspt_qru = sp_mat(locations.t(), vValues, iNs*iNp*iNt,iNq*iNr*iNu);
        Nspt_qru = 1;
        return Vspt_qru;
    }
    else{
        return Vspt_qru;
    }
}

void flexmat6::update_as_spt_qur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vt*iNs*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vr*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spt_qur(){
    if(Nspt_qur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vt*iNp*iNs;
        locations.col(1) = vq + vu*iNq + vr*iNu*iNq;
        Vspt_qur = sp_mat(locations.t(), vValues, iNs*iNp*iNt,iNq*iNu*iNr);
        Nspt_qur = 1;
        return Vspt_qur;
    }
    else{
        return Vspt_qur;
    }
}

void flexmat6::update_as_spt_rqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vt*iNs*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vu*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spt_rqu(){
    if(Nspt_rqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vt*iNp*iNs;
        locations.col(1) = vr + vq*iNr + vu*iNq*iNr;
        Vspt_rqu = sp_mat(locations.t(), vValues, iNs*iNp*iNt,iNr*iNq*iNu);
        Nspt_rqu = 1;
        return Vspt_rqu;
    }
    else{
        return Vspt_rqu;
    }
}

void flexmat6::update_as_spt_ruq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vt*iNs*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vq*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spt_ruq(){
    if(Nspt_ruq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vt*iNp*iNs;
        locations.col(1) = vr + vu*iNr + vq*iNu*iNr;
        Vspt_ruq = sp_mat(locations.t(), vValues, iNs*iNp*iNt,iNr*iNu*iNq);
        Nspt_ruq = 1;
        return Vspt_ruq;
    }
    else{
        return Vspt_ruq;
    }
}

void flexmat6::update_as_spt_uqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vt*iNs*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vr*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spt_uqr(){
    if(Nspt_uqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vt*iNp*iNs;
        locations.col(1) = vu + vq*iNu + vr*iNq*iNu;
        Vspt_uqr = sp_mat(locations.t(), vValues, iNs*iNp*iNt,iNu*iNq*iNr);
        Nspt_uqr = 1;
        return Vspt_uqr;
    }
    else{
        return Vspt_uqr;
    }
}

void flexmat6::update_as_spt_urq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vt*iNs*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vq*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spt_urq(){
    if(Nspt_urq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vt*iNp*iNs;
        locations.col(1) = vu + vr*iNu + vq*iNr*iNu;
        Vspt_urq = sp_mat(locations.t(), vValues, iNs*iNp*iNt,iNu*iNr*iNq);
        Nspt_urq = 1;
        return Vspt_urq;
    }
    else{
        return Vspt_urq;
    }
}

void flexmat6::update_as_spu_qrt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vu*iNs*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vt*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spu_qrt(){
    if(Nspu_qrt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vu*iNp*iNs;
        locations.col(1) = vq + vr*iNq + vt*iNr*iNq;
        Vspu_qrt = sp_mat(locations.t(), vValues, iNs*iNp*iNu,iNq*iNr*iNt);
        Nspu_qrt = 1;
        return Vspu_qrt;
    }
    else{
        return Vspu_qrt;
    }
}

void flexmat6::update_as_spu_qtr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vu*iNs*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vr*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spu_qtr(){
    if(Nspu_qtr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vu*iNp*iNs;
        locations.col(1) = vq + vt*iNq + vr*iNt*iNq;
        Vspu_qtr = sp_mat(locations.t(), vValues, iNs*iNp*iNu,iNq*iNt*iNr);
        Nspu_qtr = 1;
        return Vspu_qtr;
    }
    else{
        return Vspu_qtr;
    }
}

void flexmat6::update_as_spu_rqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vu*iNs*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vt*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spu_rqt(){
    if(Nspu_rqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vu*iNp*iNs;
        locations.col(1) = vr + vq*iNr + vt*iNq*iNr;
        Vspu_rqt = sp_mat(locations.t(), vValues, iNs*iNp*iNu,iNr*iNq*iNt);
        Nspu_rqt = 1;
        return Vspu_rqt;
    }
    else{
        return Vspu_rqt;
    }
}

void flexmat6::update_as_spu_rtq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vu*iNs*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vq*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spu_rtq(){
    if(Nspu_rtq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vu*iNp*iNs;
        locations.col(1) = vr + vt*iNr + vq*iNt*iNr;
        Vspu_rtq = sp_mat(locations.t(), vValues, iNs*iNp*iNu,iNr*iNt*iNq);
        Nspu_rtq = 1;
        return Vspu_rtq;
    }
    else{
        return Vspu_rtq;
    }
}

void flexmat6::update_as_spu_tqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vu*iNs*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vr*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spu_tqr(){
    if(Nspu_tqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vu*iNp*iNs;
        locations.col(1) = vt + vq*iNt + vr*iNq*iNt;
        Vspu_tqr = sp_mat(locations.t(), vValues, iNs*iNp*iNu,iNt*iNq*iNr);
        Nspu_tqr = 1;
        return Vspu_tqr;
    }
    else{
        return Vspu_tqr;
    }
}

void flexmat6::update_as_spu_trq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNs - vu*iNs*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vq*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::spu_trq(){
    if(Nspu_trq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vp*iNs + vu*iNp*iNs;
        locations.col(1) = vt + vr*iNt + vq*iNr*iNt;
        Vspu_trq = sp_mat(locations.t(), vValues, iNs*iNp*iNu,iNt*iNr*iNq);
        Nspu_trq = 1;
        return Vspu_trq;
    }
    else{
        return Vspu_trq;
    }
}

void flexmat6::update_as_sqp_rtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vp*iNs*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vu*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqp_rtu(){
    if(Nsqp_rtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vp*iNq*iNs;
        locations.col(1) = vr + vt*iNr + vu*iNt*iNr;
        Vsqp_rtu = sp_mat(locations.t(), vValues, iNs*iNq*iNp,iNr*iNt*iNu);
        Nsqp_rtu = 1;
        return Vsqp_rtu;
    }
    else{
        return Vsqp_rtu;
    }
}

void flexmat6::update_as_sqp_rut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vp*iNs*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vt*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqp_rut(){
    if(Nsqp_rut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vp*iNq*iNs;
        locations.col(1) = vr + vu*iNr + vt*iNu*iNr;
        Vsqp_rut = sp_mat(locations.t(), vValues, iNs*iNq*iNp,iNr*iNu*iNt);
        Nsqp_rut = 1;
        return Vsqp_rut;
    }
    else{
        return Vsqp_rut;
    }
}

void flexmat6::update_as_sqp_tru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vp*iNs*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vu*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqp_tru(){
    if(Nsqp_tru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vp*iNq*iNs;
        locations.col(1) = vt + vr*iNt + vu*iNr*iNt;
        Vsqp_tru = sp_mat(locations.t(), vValues, iNs*iNq*iNp,iNt*iNr*iNu);
        Nsqp_tru = 1;
        return Vsqp_tru;
    }
    else{
        return Vsqp_tru;
    }
}

void flexmat6::update_as_sqp_tur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vp*iNs*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vr*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqp_tur(){
    if(Nsqp_tur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vp*iNq*iNs;
        locations.col(1) = vt + vu*iNt + vr*iNu*iNt;
        Vsqp_tur = sp_mat(locations.t(), vValues, iNs*iNq*iNp,iNt*iNu*iNr);
        Nsqp_tur = 1;
        return Vsqp_tur;
    }
    else{
        return Vsqp_tur;
    }
}

void flexmat6::update_as_sqp_urt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vp*iNs*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vt*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqp_urt(){
    if(Nsqp_urt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vp*iNq*iNs;
        locations.col(1) = vu + vr*iNu + vt*iNr*iNu;
        Vsqp_urt = sp_mat(locations.t(), vValues, iNs*iNq*iNp,iNu*iNr*iNt);
        Nsqp_urt = 1;
        return Vsqp_urt;
    }
    else{
        return Vsqp_urt;
    }
}

void flexmat6::update_as_sqp_utr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vp*iNs*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vr*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqp_utr(){
    if(Nsqp_utr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vp*iNq*iNs;
        locations.col(1) = vu + vt*iNu + vr*iNt*iNu;
        Vsqp_utr = sp_mat(locations.t(), vValues, iNs*iNq*iNp,iNu*iNt*iNr);
        Nsqp_utr = 1;
        return Vsqp_utr;
    }
    else{
        return Vsqp_utr;
    }
}

void flexmat6::update_as_sqr_ptu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vr*iNs*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vu*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqr_ptu(){
    if(Nsqr_ptu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vr*iNq*iNs;
        locations.col(1) = vp + vt*iNp + vu*iNt*iNp;
        Vsqr_ptu = sp_mat(locations.t(), vValues, iNs*iNq*iNr,iNp*iNt*iNu);
        Nsqr_ptu = 1;
        return Vsqr_ptu;
    }
    else{
        return Vsqr_ptu;
    }
}

void flexmat6::update_as_sqr_put(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vr*iNs*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vt*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqr_put(){
    if(Nsqr_put == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vr*iNq*iNs;
        locations.col(1) = vp + vu*iNp + vt*iNu*iNp;
        Vsqr_put = sp_mat(locations.t(), vValues, iNs*iNq*iNr,iNp*iNu*iNt);
        Nsqr_put = 1;
        return Vsqr_put;
    }
    else{
        return Vsqr_put;
    }
}

void flexmat6::update_as_sqr_tpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vr*iNs*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vu*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqr_tpu(){
    if(Nsqr_tpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vr*iNq*iNs;
        locations.col(1) = vt + vp*iNt + vu*iNp*iNt;
        Vsqr_tpu = sp_mat(locations.t(), vValues, iNs*iNq*iNr,iNt*iNp*iNu);
        Nsqr_tpu = 1;
        return Vsqr_tpu;
    }
    else{
        return Vsqr_tpu;
    }
}

void flexmat6::update_as_sqr_tup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vr*iNs*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vp*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqr_tup(){
    if(Nsqr_tup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vr*iNq*iNs;
        locations.col(1) = vt + vu*iNt + vp*iNu*iNt;
        Vsqr_tup = sp_mat(locations.t(), vValues, iNs*iNq*iNr,iNt*iNu*iNp);
        Nsqr_tup = 1;
        return Vsqr_tup;
    }
    else{
        return Vsqr_tup;
    }
}

void flexmat6::update_as_sqr_upt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vr*iNs*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vt*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqr_upt(){
    if(Nsqr_upt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vr*iNq*iNs;
        locations.col(1) = vu + vp*iNu + vt*iNp*iNu;
        Vsqr_upt = sp_mat(locations.t(), vValues, iNs*iNq*iNr,iNu*iNp*iNt);
        Nsqr_upt = 1;
        return Vsqr_upt;
    }
    else{
        return Vsqr_upt;
    }
}

void flexmat6::update_as_sqr_utp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vr*iNs*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vp*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqr_utp(){
    if(Nsqr_utp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vr*iNq*iNs;
        locations.col(1) = vu + vt*iNu + vp*iNt*iNu;
        Vsqr_utp = sp_mat(locations.t(), vValues, iNs*iNq*iNr,iNu*iNt*iNp);
        Nsqr_utp = 1;
        return Vsqr_utp;
    }
    else{
        return Vsqr_utp;
    }
}

void flexmat6::update_as_sqt_pru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vt*iNs*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vu*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqt_pru(){
    if(Nsqt_pru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vt*iNq*iNs;
        locations.col(1) = vp + vr*iNp + vu*iNr*iNp;
        Vsqt_pru = sp_mat(locations.t(), vValues, iNs*iNq*iNt,iNp*iNr*iNu);
        Nsqt_pru = 1;
        return Vsqt_pru;
    }
    else{
        return Vsqt_pru;
    }
}

void flexmat6::update_as_sqt_pur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vt*iNs*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vr*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqt_pur(){
    if(Nsqt_pur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vt*iNq*iNs;
        locations.col(1) = vp + vu*iNp + vr*iNu*iNp;
        Vsqt_pur = sp_mat(locations.t(), vValues, iNs*iNq*iNt,iNp*iNu*iNr);
        Nsqt_pur = 1;
        return Vsqt_pur;
    }
    else{
        return Vsqt_pur;
    }
}

void flexmat6::update_as_sqt_rpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vt*iNs*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vu*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqt_rpu(){
    if(Nsqt_rpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vt*iNq*iNs;
        locations.col(1) = vr + vp*iNr + vu*iNp*iNr;
        Vsqt_rpu = sp_mat(locations.t(), vValues, iNs*iNq*iNt,iNr*iNp*iNu);
        Nsqt_rpu = 1;
        return Vsqt_rpu;
    }
    else{
        return Vsqt_rpu;
    }
}

void flexmat6::update_as_sqt_rup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vt*iNs*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vp*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqt_rup(){
    if(Nsqt_rup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vt*iNq*iNs;
        locations.col(1) = vr + vu*iNr + vp*iNu*iNr;
        Vsqt_rup = sp_mat(locations.t(), vValues, iNs*iNq*iNt,iNr*iNu*iNp);
        Nsqt_rup = 1;
        return Vsqt_rup;
    }
    else{
        return Vsqt_rup;
    }
}

void flexmat6::update_as_sqt_upr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vt*iNs*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vr*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqt_upr(){
    if(Nsqt_upr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vt*iNq*iNs;
        locations.col(1) = vu + vp*iNu + vr*iNp*iNu;
        Vsqt_upr = sp_mat(locations.t(), vValues, iNs*iNq*iNt,iNu*iNp*iNr);
        Nsqt_upr = 1;
        return Vsqt_upr;
    }
    else{
        return Vsqt_upr;
    }
}

void flexmat6::update_as_sqt_urp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vt*iNs*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vp*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sqt_urp(){
    if(Nsqt_urp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vt*iNq*iNs;
        locations.col(1) = vu + vr*iNu + vp*iNr*iNu;
        Vsqt_urp = sp_mat(locations.t(), vValues, iNs*iNq*iNt,iNu*iNr*iNp);
        Nsqt_urp = 1;
        return Vsqt_urp;
    }
    else{
        return Vsqt_urp;
    }
}

void flexmat6::update_as_squ_prt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vu*iNs*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vt*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::squ_prt(){
    if(Nsqu_prt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vu*iNq*iNs;
        locations.col(1) = vp + vr*iNp + vt*iNr*iNp;
        Vsqu_prt = sp_mat(locations.t(), vValues, iNs*iNq*iNu,iNp*iNr*iNt);
        Nsqu_prt = 1;
        return Vsqu_prt;
    }
    else{
        return Vsqu_prt;
    }
}

void flexmat6::update_as_squ_ptr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vu*iNs*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vr*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::squ_ptr(){
    if(Nsqu_ptr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vu*iNq*iNs;
        locations.col(1) = vp + vt*iNp + vr*iNt*iNp;
        Vsqu_ptr = sp_mat(locations.t(), vValues, iNs*iNq*iNu,iNp*iNt*iNr);
        Nsqu_ptr = 1;
        return Vsqu_ptr;
    }
    else{
        return Vsqu_ptr;
    }
}

void flexmat6::update_as_squ_rpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vu*iNs*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vt*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::squ_rpt(){
    if(Nsqu_rpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vu*iNq*iNs;
        locations.col(1) = vr + vp*iNr + vt*iNp*iNr;
        Vsqu_rpt = sp_mat(locations.t(), vValues, iNs*iNq*iNu,iNr*iNp*iNt);
        Nsqu_rpt = 1;
        return Vsqu_rpt;
    }
    else{
        return Vsqu_rpt;
    }
}

void flexmat6::update_as_squ_rtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vu*iNs*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vp*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::squ_rtp(){
    if(Nsqu_rtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vu*iNq*iNs;
        locations.col(1) = vr + vt*iNr + vp*iNt*iNr;
        Vsqu_rtp = sp_mat(locations.t(), vValues, iNs*iNq*iNu,iNr*iNt*iNp);
        Nsqu_rtp = 1;
        return Vsqu_rtp;
    }
    else{
        return Vsqu_rtp;
    }
}

void flexmat6::update_as_squ_tpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vu*iNs*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vr*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::squ_tpr(){
    if(Nsqu_tpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vu*iNq*iNs;
        locations.col(1) = vt + vp*iNt + vr*iNp*iNt;
        Vsqu_tpr = sp_mat(locations.t(), vValues, iNs*iNq*iNu,iNt*iNp*iNr);
        Nsqu_tpr = 1;
        return Vsqu_tpr;
    }
    else{
        return Vsqu_tpr;
    }
}

void flexmat6::update_as_squ_trp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNs - vu*iNs*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vp*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::squ_trp(){
    if(Nsqu_trp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vq*iNs + vu*iNq*iNs;
        locations.col(1) = vt + vr*iNt + vp*iNr*iNt;
        Vsqu_trp = sp_mat(locations.t(), vValues, iNs*iNq*iNu,iNt*iNr*iNp);
        Nsqu_trp = 1;
        return Vsqu_trp;
    }
    else{
        return Vsqu_trp;
    }
}

void flexmat6::update_as_srp_qtu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vp*iNs*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vu*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srp_qtu(){
    if(Nsrp_qtu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vp*iNr*iNs;
        locations.col(1) = vq + vt*iNq + vu*iNt*iNq;
        Vsrp_qtu = sp_mat(locations.t(), vValues, iNs*iNr*iNp,iNq*iNt*iNu);
        Nsrp_qtu = 1;
        return Vsrp_qtu;
    }
    else{
        return Vsrp_qtu;
    }
}

void flexmat6::update_as_srp_qut(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vp*iNs*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vt*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srp_qut(){
    if(Nsrp_qut == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vp*iNr*iNs;
        locations.col(1) = vq + vu*iNq + vt*iNu*iNq;
        Vsrp_qut = sp_mat(locations.t(), vValues, iNs*iNr*iNp,iNq*iNu*iNt);
        Nsrp_qut = 1;
        return Vsrp_qut;
    }
    else{
        return Vsrp_qut;
    }
}

void flexmat6::update_as_srp_tqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vp*iNs*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vu*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srp_tqu(){
    if(Nsrp_tqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vp*iNr*iNs;
        locations.col(1) = vt + vq*iNt + vu*iNq*iNt;
        Vsrp_tqu = sp_mat(locations.t(), vValues, iNs*iNr*iNp,iNt*iNq*iNu);
        Nsrp_tqu = 1;
        return Vsrp_tqu;
    }
    else{
        return Vsrp_tqu;
    }
}

void flexmat6::update_as_srp_tuq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vp*iNs*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vq*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srp_tuq(){
    if(Nsrp_tuq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vp*iNr*iNs;
        locations.col(1) = vt + vu*iNt + vq*iNu*iNt;
        Vsrp_tuq = sp_mat(locations.t(), vValues, iNs*iNr*iNp,iNt*iNu*iNq);
        Nsrp_tuq = 1;
        return Vsrp_tuq;
    }
    else{
        return Vsrp_tuq;
    }
}

void flexmat6::update_as_srp_uqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vp*iNs*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vt*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srp_uqt(){
    if(Nsrp_uqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vp*iNr*iNs;
        locations.col(1) = vu + vq*iNu + vt*iNq*iNu;
        Vsrp_uqt = sp_mat(locations.t(), vValues, iNs*iNr*iNp,iNu*iNq*iNt);
        Nsrp_uqt = 1;
        return Vsrp_uqt;
    }
    else{
        return Vsrp_uqt;
    }
}

void flexmat6::update_as_srp_utq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vp*iNs*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vq*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srp_utq(){
    if(Nsrp_utq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vp*iNr*iNs;
        locations.col(1) = vu + vt*iNu + vq*iNt*iNu;
        Vsrp_utq = sp_mat(locations.t(), vValues, iNs*iNr*iNp,iNu*iNt*iNq);
        Nsrp_utq = 1;
        return Vsrp_utq;
    }
    else{
        return Vsrp_utq;
    }
}

void flexmat6::update_as_srq_ptu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vq*iNs*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vu*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srq_ptu(){
    if(Nsrq_ptu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vq*iNr*iNs;
        locations.col(1) = vp + vt*iNp + vu*iNt*iNp;
        Vsrq_ptu = sp_mat(locations.t(), vValues, iNs*iNr*iNq,iNp*iNt*iNu);
        Nsrq_ptu = 1;
        return Vsrq_ptu;
    }
    else{
        return Vsrq_ptu;
    }
}

void flexmat6::update_as_srq_put(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vq*iNs*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vt*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srq_put(){
    if(Nsrq_put == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vq*iNr*iNs;
        locations.col(1) = vp + vu*iNp + vt*iNu*iNp;
        Vsrq_put = sp_mat(locations.t(), vValues, iNs*iNr*iNq,iNp*iNu*iNt);
        Nsrq_put = 1;
        return Vsrq_put;
    }
    else{
        return Vsrq_put;
    }
}

void flexmat6::update_as_srq_tpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vq*iNs*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vu*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srq_tpu(){
    if(Nsrq_tpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vq*iNr*iNs;
        locations.col(1) = vt + vp*iNt + vu*iNp*iNt;
        Vsrq_tpu = sp_mat(locations.t(), vValues, iNs*iNr*iNq,iNt*iNp*iNu);
        Nsrq_tpu = 1;
        return Vsrq_tpu;
    }
    else{
        return Vsrq_tpu;
    }
}

void flexmat6::update_as_srq_tup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vq*iNs*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vu*iNt - vp*iNt*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srq_tup(){
    if(Nsrq_tup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vq*iNr*iNs;
        locations.col(1) = vt + vu*iNt + vp*iNu*iNt;
        Vsrq_tup = sp_mat(locations.t(), vValues, iNs*iNr*iNq,iNt*iNu*iNp);
        Nsrq_tup = 1;
        return Vsrq_tup;
    }
    else{
        return Vsrq_tup;
    }
}

void flexmat6::update_as_srq_upt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vq*iNs*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vt*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srq_upt(){
    if(Nsrq_upt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vq*iNr*iNs;
        locations.col(1) = vu + vp*iNu + vt*iNp*iNu;
        Vsrq_upt = sp_mat(locations.t(), vValues, iNs*iNr*iNq,iNu*iNp*iNt);
        Nsrq_upt = 1;
        return Vsrq_upt;
    }
    else{
        return Vsrq_upt;
    }
}

void flexmat6::update_as_srq_utp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vq*iNs*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vt*iNu - vp*iNu*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srq_utp(){
    if(Nsrq_utp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vq*iNr*iNs;
        locations.col(1) = vu + vt*iNu + vp*iNt*iNu;
        Vsrq_utp = sp_mat(locations.t(), vValues, iNs*iNr*iNq,iNu*iNt*iNp);
        Nsrq_utp = 1;
        return Vsrq_utp;
    }
    else{
        return Vsrq_utp;
    }
}

void flexmat6::update_as_srt_pqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vt*iNs*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vu*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srt_pqu(){
    if(Nsrt_pqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vt*iNr*iNs;
        locations.col(1) = vp + vq*iNp + vu*iNq*iNp;
        Vsrt_pqu = sp_mat(locations.t(), vValues, iNs*iNr*iNt,iNp*iNq*iNu);
        Nsrt_pqu = 1;
        return Vsrt_pqu;
    }
    else{
        return Vsrt_pqu;
    }
}

void flexmat6::update_as_srt_puq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vt*iNs*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vq*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srt_puq(){
    if(Nsrt_puq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vt*iNr*iNs;
        locations.col(1) = vp + vu*iNp + vq*iNu*iNp;
        Vsrt_puq = sp_mat(locations.t(), vValues, iNs*iNr*iNt,iNp*iNu*iNq);
        Nsrt_puq = 1;
        return Vsrt_puq;
    }
    else{
        return Vsrt_puq;
    }
}

void flexmat6::update_as_srt_qpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vt*iNs*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vu*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srt_qpu(){
    if(Nsrt_qpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vt*iNr*iNs;
        locations.col(1) = vq + vp*iNq + vu*iNp*iNq;
        Vsrt_qpu = sp_mat(locations.t(), vValues, iNs*iNr*iNt,iNq*iNp*iNu);
        Nsrt_qpu = 1;
        return Vsrt_qpu;
    }
    else{
        return Vsrt_qpu;
    }
}

void flexmat6::update_as_srt_qup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vt*iNs*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vp*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srt_qup(){
    if(Nsrt_qup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vt*iNr*iNs;
        locations.col(1) = vq + vu*iNq + vp*iNu*iNq;
        Vsrt_qup = sp_mat(locations.t(), vValues, iNs*iNr*iNt,iNq*iNu*iNp);
        Nsrt_qup = 1;
        return Vsrt_qup;
    }
    else{
        return Vsrt_qup;
    }
}

void flexmat6::update_as_srt_upq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vt*iNs*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vq*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srt_upq(){
    if(Nsrt_upq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vt*iNr*iNs;
        locations.col(1) = vu + vp*iNu + vq*iNp*iNu;
        Vsrt_upq = sp_mat(locations.t(), vValues, iNs*iNr*iNt,iNu*iNp*iNq);
        Nsrt_upq = 1;
        return Vsrt_upq;
    }
    else{
        return Vsrt_upq;
    }
}

void flexmat6::update_as_srt_uqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vt*iNs*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vp*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::srt_uqp(){
    if(Nsrt_uqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vt*iNr*iNs;
        locations.col(1) = vu + vq*iNu + vp*iNq*iNu;
        Vsrt_uqp = sp_mat(locations.t(), vValues, iNs*iNr*iNt,iNu*iNq*iNp);
        Nsrt_uqp = 1;
        return Vsrt_uqp;
    }
    else{
        return Vsrt_uqp;
    }
}

void flexmat6::update_as_sru_pqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vu*iNs*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vt*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sru_pqt(){
    if(Nsru_pqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vu*iNr*iNs;
        locations.col(1) = vp + vq*iNp + vt*iNq*iNp;
        Vsru_pqt = sp_mat(locations.t(), vValues, iNs*iNr*iNu,iNp*iNq*iNt);
        Nsru_pqt = 1;
        return Vsru_pqt;
    }
    else{
        return Vsru_pqt;
    }
}

void flexmat6::update_as_sru_ptq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vu*iNs*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vq*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sru_ptq(){
    if(Nsru_ptq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vu*iNr*iNs;
        locations.col(1) = vp + vt*iNp + vq*iNt*iNp;
        Vsru_ptq = sp_mat(locations.t(), vValues, iNs*iNr*iNu,iNp*iNt*iNq);
        Nsru_ptq = 1;
        return Vsru_ptq;
    }
    else{
        return Vsru_ptq;
    }
}

void flexmat6::update_as_sru_qpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vu*iNs*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vt*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sru_qpt(){
    if(Nsru_qpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vu*iNr*iNs;
        locations.col(1) = vq + vp*iNq + vt*iNp*iNq;
        Vsru_qpt = sp_mat(locations.t(), vValues, iNs*iNr*iNu,iNq*iNp*iNt);
        Nsru_qpt = 1;
        return Vsru_qpt;
    }
    else{
        return Vsru_qpt;
    }
}

void flexmat6::update_as_sru_qtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vu*iNs*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vp*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sru_qtp(){
    if(Nsru_qtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vu*iNr*iNs;
        locations.col(1) = vq + vt*iNq + vp*iNt*iNq;
        Vsru_qtp = sp_mat(locations.t(), vValues, iNs*iNr*iNu,iNq*iNt*iNp);
        Nsru_qtp = 1;
        return Vsru_qtp;
    }
    else{
        return Vsru_qtp;
    }
}

void flexmat6::update_as_sru_tpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vu*iNs*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vq*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sru_tpq(){
    if(Nsru_tpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vu*iNr*iNs;
        locations.col(1) = vt + vp*iNt + vq*iNp*iNt;
        Vsru_tpq = sp_mat(locations.t(), vValues, iNs*iNr*iNu,iNt*iNp*iNq);
        Nsru_tpq = 1;
        return Vsru_tpq;
    }
    else{
        return Vsru_tpq;
    }
}

void flexmat6::update_as_sru_tqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNs - vu*iNs*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vp*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sru_tqp(){
    if(Nsru_tqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vr*iNs + vu*iNr*iNs;
        locations.col(1) = vt + vq*iNt + vp*iNq*iNt;
        Vsru_tqp = sp_mat(locations.t(), vValues, iNs*iNr*iNu,iNt*iNq*iNp);
        Nsru_tqp = 1;
        return Vsru_tqp;
    }
    else{
        return Vsru_tqp;
    }
}

void flexmat6::update_as_stp_qru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vp*iNs*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vu*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stp_qru(){
    if(Nstp_qru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vp*iNt*iNs;
        locations.col(1) = vq + vr*iNq + vu*iNr*iNq;
        Vstp_qru = sp_mat(locations.t(), vValues, iNs*iNt*iNp,iNq*iNr*iNu);
        Nstp_qru = 1;
        return Vstp_qru;
    }
    else{
        return Vstp_qru;
    }
}

void flexmat6::update_as_stp_qur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vp*iNs*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vr*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stp_qur(){
    if(Nstp_qur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vp*iNt*iNs;
        locations.col(1) = vq + vu*iNq + vr*iNu*iNq;
        Vstp_qur = sp_mat(locations.t(), vValues, iNs*iNt*iNp,iNq*iNu*iNr);
        Nstp_qur = 1;
        return Vstp_qur;
    }
    else{
        return Vstp_qur;
    }
}

void flexmat6::update_as_stp_rqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vp*iNs*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vu*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stp_rqu(){
    if(Nstp_rqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vp*iNt*iNs;
        locations.col(1) = vr + vq*iNr + vu*iNq*iNr;
        Vstp_rqu = sp_mat(locations.t(), vValues, iNs*iNt*iNp,iNr*iNq*iNu);
        Nstp_rqu = 1;
        return Vstp_rqu;
    }
    else{
        return Vstp_rqu;
    }
}

void flexmat6::update_as_stp_ruq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vp*iNs*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vq*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stp_ruq(){
    if(Nstp_ruq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vp*iNt*iNs;
        locations.col(1) = vr + vu*iNr + vq*iNu*iNr;
        Vstp_ruq = sp_mat(locations.t(), vValues, iNs*iNt*iNp,iNr*iNu*iNq);
        Nstp_ruq = 1;
        return Vstp_ruq;
    }
    else{
        return Vstp_ruq;
    }
}

void flexmat6::update_as_stp_uqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vp*iNs*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vr*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stp_uqr(){
    if(Nstp_uqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vp*iNt*iNs;
        locations.col(1) = vu + vq*iNu + vr*iNq*iNu;
        Vstp_uqr = sp_mat(locations.t(), vValues, iNs*iNt*iNp,iNu*iNq*iNr);
        Nstp_uqr = 1;
        return Vstp_uqr;
    }
    else{
        return Vstp_uqr;
    }
}

void flexmat6::update_as_stp_urq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vp*iNs*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vq*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stp_urq(){
    if(Nstp_urq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vp*iNt*iNs;
        locations.col(1) = vu + vr*iNu + vq*iNr*iNu;
        Vstp_urq = sp_mat(locations.t(), vValues, iNs*iNt*iNp,iNu*iNr*iNq);
        Nstp_urq = 1;
        return Vstp_urq;
    }
    else{
        return Vstp_urq;
    }
}

void flexmat6::update_as_stq_pru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vq*iNs*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vu*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stq_pru(){
    if(Nstq_pru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vq*iNt*iNs;
        locations.col(1) = vp + vr*iNp + vu*iNr*iNp;
        Vstq_pru = sp_mat(locations.t(), vValues, iNs*iNt*iNq,iNp*iNr*iNu);
        Nstq_pru = 1;
        return Vstq_pru;
    }
    else{
        return Vstq_pru;
    }
}

void flexmat6::update_as_stq_pur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vq*iNs*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vr*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stq_pur(){
    if(Nstq_pur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vq*iNt*iNs;
        locations.col(1) = vp + vu*iNp + vr*iNu*iNp;
        Vstq_pur = sp_mat(locations.t(), vValues, iNs*iNt*iNq,iNp*iNu*iNr);
        Nstq_pur = 1;
        return Vstq_pur;
    }
    else{
        return Vstq_pur;
    }
}

void flexmat6::update_as_stq_rpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vq*iNs*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vu*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stq_rpu(){
    if(Nstq_rpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vq*iNt*iNs;
        locations.col(1) = vr + vp*iNr + vu*iNp*iNr;
        Vstq_rpu = sp_mat(locations.t(), vValues, iNs*iNt*iNq,iNr*iNp*iNu);
        Nstq_rpu = 1;
        return Vstq_rpu;
    }
    else{
        return Vstq_rpu;
    }
}

void flexmat6::update_as_stq_rup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vq*iNs*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vp*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stq_rup(){
    if(Nstq_rup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vq*iNt*iNs;
        locations.col(1) = vr + vu*iNr + vp*iNu*iNr;
        Vstq_rup = sp_mat(locations.t(), vValues, iNs*iNt*iNq,iNr*iNu*iNp);
        Nstq_rup = 1;
        return Vstq_rup;
    }
    else{
        return Vstq_rup;
    }
}

void flexmat6::update_as_stq_upr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vq*iNs*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vr*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stq_upr(){
    if(Nstq_upr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vq*iNt*iNs;
        locations.col(1) = vu + vp*iNu + vr*iNp*iNu;
        Vstq_upr = sp_mat(locations.t(), vValues, iNs*iNt*iNq,iNu*iNp*iNr);
        Nstq_upr = 1;
        return Vstq_upr;
    }
    else{
        return Vstq_upr;
    }
}

void flexmat6::update_as_stq_urp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vq*iNs*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vp*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stq_urp(){
    if(Nstq_urp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vq*iNt*iNs;
        locations.col(1) = vu + vr*iNu + vp*iNr*iNu;
        Vstq_urp = sp_mat(locations.t(), vValues, iNs*iNt*iNq,iNu*iNr*iNp);
        Nstq_urp = 1;
        return Vstq_urp;
    }
    else{
        return Vstq_urp;
    }
}

void flexmat6::update_as_str_pqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vr*iNs*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vu*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::str_pqu(){
    if(Nstr_pqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vr*iNt*iNs;
        locations.col(1) = vp + vq*iNp + vu*iNq*iNp;
        Vstr_pqu = sp_mat(locations.t(), vValues, iNs*iNt*iNr,iNp*iNq*iNu);
        Nstr_pqu = 1;
        return Vstr_pqu;
    }
    else{
        return Vstr_pqu;
    }
}

void flexmat6::update_as_str_puq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vr*iNs*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vq*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::str_puq(){
    if(Nstr_puq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vr*iNt*iNs;
        locations.col(1) = vp + vu*iNp + vq*iNu*iNp;
        Vstr_puq = sp_mat(locations.t(), vValues, iNs*iNt*iNr,iNp*iNu*iNq);
        Nstr_puq = 1;
        return Vstr_puq;
    }
    else{
        return Vstr_puq;
    }
}

void flexmat6::update_as_str_qpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vr*iNs*iNt)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vu*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::str_qpu(){
    if(Nstr_qpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vr*iNt*iNs;
        locations.col(1) = vq + vp*iNq + vu*iNp*iNq;
        Vstr_qpu = sp_mat(locations.t(), vValues, iNs*iNt*iNr,iNq*iNp*iNu);
        Nstr_qpu = 1;
        return Vstr_qpu;
    }
    else{
        return Vstr_qpu;
    }
}

void flexmat6::update_as_str_qup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vr*iNs*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vp*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::str_qup(){
    if(Nstr_qup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vr*iNt*iNs;
        locations.col(1) = vq + vu*iNq + vp*iNu*iNq;
        Vstr_qup = sp_mat(locations.t(), vValues, iNs*iNt*iNr,iNq*iNu*iNp);
        Nstr_qup = 1;
        return Vstr_qup;
    }
    else{
        return Vstr_qup;
    }
}

void flexmat6::update_as_str_upq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vr*iNs*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vq*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::str_upq(){
    if(Nstr_upq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vr*iNt*iNs;
        locations.col(1) = vu + vp*iNu + vq*iNp*iNu;
        Vstr_upq = sp_mat(locations.t(), vValues, iNs*iNt*iNr,iNu*iNp*iNq);
        Nstr_upq = 1;
        return Vstr_upq;
    }
    else{
        return Vstr_upq;
    }
}

void flexmat6::update_as_str_uqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vr*iNs*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vp*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::str_uqp(){
    if(Nstr_uqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vr*iNt*iNs;
        locations.col(1) = vu + vq*iNu + vp*iNq*iNu;
        Vstr_uqp = sp_mat(locations.t(), vValues, iNs*iNt*iNr,iNu*iNq*iNp);
        Nstr_uqp = 1;
        return Vstr_uqp;
    }
    else{
        return Vstr_uqp;
    }
}

void flexmat6::update_as_stu_pqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vu*iNs*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vr*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stu_pqr(){
    if(Nstu_pqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vu*iNt*iNs;
        locations.col(1) = vp + vq*iNp + vr*iNq*iNp;
        Vstu_pqr = sp_mat(locations.t(), vValues, iNs*iNt*iNu,iNp*iNq*iNr);
        Nstu_pqr = 1;
        return Vstu_pqr;
    }
    else{
        return Vstu_pqr;
    }
}

void flexmat6::update_as_stu_prq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vu*iNs*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vq*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stu_prq(){
    if(Nstu_prq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vu*iNt*iNs;
        locations.col(1) = vp + vr*iNp + vq*iNr*iNp;
        Vstu_prq = sp_mat(locations.t(), vValues, iNs*iNt*iNu,iNp*iNr*iNq);
        Nstu_prq = 1;
        return Vstu_prq;
    }
    else{
        return Vstu_prq;
    }
}

void flexmat6::update_as_stu_qpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vu*iNs*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vr*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stu_qpr(){
    if(Nstu_qpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vu*iNt*iNs;
        locations.col(1) = vq + vp*iNq + vr*iNp*iNq;
        Vstu_qpr = sp_mat(locations.t(), vValues, iNs*iNt*iNu,iNq*iNp*iNr);
        Nstu_qpr = 1;
        return Vstu_qpr;
    }
    else{
        return Vstu_qpr;
    }
}

void flexmat6::update_as_stu_qrp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vu*iNs*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vp*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stu_qrp(){
    if(Nstu_qrp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vu*iNt*iNs;
        locations.col(1) = vq + vr*iNq + vp*iNr*iNq;
        Vstu_qrp = sp_mat(locations.t(), vValues, iNs*iNt*iNu,iNq*iNr*iNp);
        Nstu_qrp = 1;
        return Vstu_qrp;
    }
    else{
        return Vstu_qrp;
    }
}

void flexmat6::update_as_stu_rpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vu*iNs*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vq*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stu_rpq(){
    if(Nstu_rpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vu*iNt*iNs;
        locations.col(1) = vr + vp*iNr + vq*iNp*iNr;
        Vstu_rpq = sp_mat(locations.t(), vValues, iNs*iNt*iNu,iNr*iNp*iNq);
        Nstu_rpq = 1;
        return Vstu_rpq;
    }
    else{
        return Vstu_rpq;
    }
}

void flexmat6::update_as_stu_rqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNs - vu*iNs*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vp*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::stu_rqp(){
    if(Nstu_rqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vt*iNs + vu*iNt*iNs;
        locations.col(1) = vr + vq*iNr + vp*iNq*iNr;
        Vstu_rqp = sp_mat(locations.t(), vValues, iNs*iNt*iNu,iNr*iNq*iNp);
        Nstu_rqp = 1;
        return Vstu_rqp;
    }
    else{
        return Vstu_rqp;
    }
}

void flexmat6::update_as_sup_qrt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vp*iNs*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vt*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sup_qrt(){
    if(Nsup_qrt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vp*iNu*iNs;
        locations.col(1) = vq + vr*iNq + vt*iNr*iNq;
        Vsup_qrt = sp_mat(locations.t(), vValues, iNs*iNu*iNp,iNq*iNr*iNt);
        Nsup_qrt = 1;
        return Vsup_qrt;
    }
    else{
        return Vsup_qrt;
    }
}

void flexmat6::update_as_sup_qtr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vp*iNs*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vr*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sup_qtr(){
    if(Nsup_qtr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vp*iNu*iNs;
        locations.col(1) = vq + vt*iNq + vr*iNt*iNq;
        Vsup_qtr = sp_mat(locations.t(), vValues, iNs*iNu*iNp,iNq*iNt*iNr);
        Nsup_qtr = 1;
        return Vsup_qtr;
    }
    else{
        return Vsup_qtr;
    }
}

void flexmat6::update_as_sup_rqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vp*iNs*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vt*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sup_rqt(){
    if(Nsup_rqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vp*iNu*iNs;
        locations.col(1) = vr + vq*iNr + vt*iNq*iNr;
        Vsup_rqt = sp_mat(locations.t(), vValues, iNs*iNu*iNp,iNr*iNq*iNt);
        Nsup_rqt = 1;
        return Vsup_rqt;
    }
    else{
        return Vsup_rqt;
    }
}

void flexmat6::update_as_sup_rtq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vp*iNs*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vq*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sup_rtq(){
    if(Nsup_rtq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vp*iNu*iNs;
        locations.col(1) = vr + vt*iNr + vq*iNt*iNr;
        Vsup_rtq = sp_mat(locations.t(), vValues, iNs*iNu*iNp,iNr*iNt*iNq);
        Nsup_rtq = 1;
        return Vsup_rtq;
    }
    else{
        return Vsup_rtq;
    }
}

void flexmat6::update_as_sup_tqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vp*iNs*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vr*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sup_tqr(){
    if(Nsup_tqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vp*iNu*iNs;
        locations.col(1) = vt + vq*iNt + vr*iNq*iNt;
        Vsup_tqr = sp_mat(locations.t(), vValues, iNs*iNu*iNp,iNt*iNq*iNr);
        Nsup_tqr = 1;
        return Vsup_tqr;
    }
    else{
        return Vsup_tqr;
    }
}

void flexmat6::update_as_sup_trq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vp*iNs*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vq*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sup_trq(){
    if(Nsup_trq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vp*iNu*iNs;
        locations.col(1) = vt + vr*iNt + vq*iNr*iNt;
        Vsup_trq = sp_mat(locations.t(), vValues, iNs*iNu*iNp,iNt*iNr*iNq);
        Nsup_trq = 1;
        return Vsup_trq;
    }
    else{
        return Vsup_trq;
    }
}

void flexmat6::update_as_suq_prt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vq*iNs*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vt*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::suq_prt(){
    if(Nsuq_prt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vq*iNu*iNs;
        locations.col(1) = vp + vr*iNp + vt*iNr*iNp;
        Vsuq_prt = sp_mat(locations.t(), vValues, iNs*iNu*iNq,iNp*iNr*iNt);
        Nsuq_prt = 1;
        return Vsuq_prt;
    }
    else{
        return Vsuq_prt;
    }
}

void flexmat6::update_as_suq_ptr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vq*iNs*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vr*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::suq_ptr(){
    if(Nsuq_ptr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vq*iNu*iNs;
        locations.col(1) = vp + vt*iNp + vr*iNt*iNp;
        Vsuq_ptr = sp_mat(locations.t(), vValues, iNs*iNu*iNq,iNp*iNt*iNr);
        Nsuq_ptr = 1;
        return Vsuq_ptr;
    }
    else{
        return Vsuq_ptr;
    }
}

void flexmat6::update_as_suq_rpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vq*iNs*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vt*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::suq_rpt(){
    if(Nsuq_rpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vq*iNu*iNs;
        locations.col(1) = vr + vp*iNr + vt*iNp*iNr;
        Vsuq_rpt = sp_mat(locations.t(), vValues, iNs*iNu*iNq,iNr*iNp*iNt);
        Nsuq_rpt = 1;
        return Vsuq_rpt;
    }
    else{
        return Vsuq_rpt;
    }
}

void flexmat6::update_as_suq_rtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vq*iNs*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vp*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::suq_rtp(){
    if(Nsuq_rtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vq*iNu*iNs;
        locations.col(1) = vr + vt*iNr + vp*iNt*iNr;
        Vsuq_rtp = sp_mat(locations.t(), vValues, iNs*iNu*iNq,iNr*iNt*iNp);
        Nsuq_rtp = 1;
        return Vsuq_rtp;
    }
    else{
        return Vsuq_rtp;
    }
}

void flexmat6::update_as_suq_tpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vq*iNs*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vr*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::suq_tpr(){
    if(Nsuq_tpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vq*iNu*iNs;
        locations.col(1) = vt + vp*iNt + vr*iNp*iNt;
        Vsuq_tpr = sp_mat(locations.t(), vValues, iNs*iNu*iNq,iNt*iNp*iNr);
        Nsuq_tpr = 1;
        return Vsuq_tpr;
    }
    else{
        return Vsuq_tpr;
    }
}

void flexmat6::update_as_suq_trp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vq*iNs*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vp*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::suq_trp(){
    if(Nsuq_trp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vq*iNu*iNs;
        locations.col(1) = vt + vr*iNt + vp*iNr*iNt;
        Vsuq_trp = sp_mat(locations.t(), vValues, iNs*iNu*iNq,iNt*iNr*iNp);
        Nsuq_trp = 1;
        return Vsuq_trp;
    }
    else{
        return Vsuq_trp;
    }
}

void flexmat6::update_as_sur_pqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vr*iNs*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vt*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sur_pqt(){
    if(Nsur_pqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vr*iNu*iNs;
        locations.col(1) = vp + vq*iNp + vt*iNq*iNp;
        Vsur_pqt = sp_mat(locations.t(), vValues, iNs*iNu*iNr,iNp*iNq*iNt);
        Nsur_pqt = 1;
        return Vsur_pqt;
    }
    else{
        return Vsur_pqt;
    }
}

void flexmat6::update_as_sur_ptq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vr*iNs*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vq*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sur_ptq(){
    if(Nsur_ptq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vr*iNu*iNs;
        locations.col(1) = vp + vt*iNp + vq*iNt*iNp;
        Vsur_ptq = sp_mat(locations.t(), vValues, iNs*iNu*iNr,iNp*iNt*iNq);
        Nsur_ptq = 1;
        return Vsur_ptq;
    }
    else{
        return Vsur_ptq;
    }
}

void flexmat6::update_as_sur_qpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vr*iNs*iNu)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vt*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sur_qpt(){
    if(Nsur_qpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vr*iNu*iNs;
        locations.col(1) = vq + vp*iNq + vt*iNp*iNq;
        Vsur_qpt = sp_mat(locations.t(), vValues, iNs*iNu*iNr,iNq*iNp*iNt);
        Nsur_qpt = 1;
        return Vsur_qpt;
    }
    else{
        return Vsur_qpt;
    }
}

void flexmat6::update_as_sur_qtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vr*iNs*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vp*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sur_qtp(){
    if(Nsur_qtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vr*iNu*iNs;
        locations.col(1) = vq + vt*iNq + vp*iNt*iNq;
        Vsur_qtp = sp_mat(locations.t(), vValues, iNs*iNu*iNr,iNq*iNt*iNp);
        Nsur_qtp = 1;
        return Vsur_qtp;
    }
    else{
        return Vsur_qtp;
    }
}

void flexmat6::update_as_sur_tpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vr*iNs*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vq*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sur_tpq(){
    if(Nsur_tpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vr*iNu*iNs;
        locations.col(1) = vt + vp*iNt + vq*iNp*iNt;
        Vsur_tpq = sp_mat(locations.t(), vValues, iNs*iNu*iNr,iNt*iNp*iNq);
        Nsur_tpq = 1;
        return Vsur_tpq;
    }
    else{
        return Vsur_tpq;
    }
}

void flexmat6::update_as_sur_tqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vr*iNs*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vp*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sur_tqp(){
    if(Nsur_tqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vr*iNu*iNs;
        locations.col(1) = vt + vq*iNt + vp*iNq*iNt;
        Vsur_tqp = sp_mat(locations.t(), vValues, iNs*iNu*iNr,iNt*iNq*iNp);
        Nsur_tqp = 1;
        return Vsur_tqp;
    }
    else{
        return Vsur_tqp;
    }
}

void flexmat6::update_as_sut_pqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vt*iNs*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vr*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sut_pqr(){
    if(Nsut_pqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vt*iNu*iNs;
        locations.col(1) = vp + vq*iNp + vr*iNq*iNp;
        Vsut_pqr = sp_mat(locations.t(), vValues, iNs*iNu*iNt,iNp*iNq*iNr);
        Nsut_pqr = 1;
        return Vsut_pqr;
    }
    else{
        return Vsut_pqr;
    }
}

void flexmat6::update_as_sut_prq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vt*iNs*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vq*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sut_prq(){
    if(Nsut_prq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vt*iNu*iNs;
        locations.col(1) = vp + vr*iNp + vq*iNr*iNp;
        Vsut_prq = sp_mat(locations.t(), vValues, iNs*iNu*iNt,iNp*iNr*iNq);
        Nsut_prq = 1;
        return Vsut_prq;
    }
    else{
        return Vsut_prq;
    }
}

void flexmat6::update_as_sut_qpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vt*iNs*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vr*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sut_qpr(){
    if(Nsut_qpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vt*iNu*iNs;
        locations.col(1) = vq + vp*iNq + vr*iNp*iNq;
        Vsut_qpr = sp_mat(locations.t(), vValues, iNs*iNu*iNt,iNq*iNp*iNr);
        Nsut_qpr = 1;
        return Vsut_qpr;
    }
    else{
        return Vsut_qpr;
    }
}

void flexmat6::update_as_sut_qrp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vt*iNs*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vp*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sut_qrp(){
    if(Nsut_qrp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vt*iNu*iNs;
        locations.col(1) = vq + vr*iNq + vp*iNr*iNq;
        Vsut_qrp = sp_mat(locations.t(), vValues, iNs*iNu*iNt,iNq*iNr*iNp);
        Nsut_qrp = 1;
        return Vsut_qrp;
    }
    else{
        return Vsut_qrp;
    }
}

void flexmat6::update_as_sut_rpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vt*iNs*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vq*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sut_rpq(){
    if(Nsut_rpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vt*iNu*iNs;
        locations.col(1) = vr + vp*iNr + vq*iNp*iNr;
        Vsut_rpq = sp_mat(locations.t(), vValues, iNs*iNu*iNt,iNr*iNp*iNq);
        Nsut_rpq = 1;
        return Vsut_rpq;
    }
    else{
        return Vsut_rpq;
    }
}

void flexmat6::update_as_sut_rqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNs - vt*iNs*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vp*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::sut_rqp(){
    if(Nsut_rqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vs + vu*iNs + vt*iNu*iNs;
        locations.col(1) = vr + vq*iNr + vp*iNq*iNr;
        Vsut_rqp = sp_mat(locations.t(), vValues, iNs*iNu*iNt,iNr*iNq*iNp);
        Nsut_rqp = 1;
        return Vsut_rqp;
    }
    else{
        return Vsut_rqp;
    }
}

void flexmat6::update_as_tpq_rsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vq*iNt*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vu*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpq_rsu(){
    if(Ntpq_rsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vq*iNp*iNt;
        locations.col(1) = vr + vs*iNr + vu*iNs*iNr;
        Vtpq_rsu = sp_mat(locations.t(), vValues, iNt*iNp*iNq,iNr*iNs*iNu);
        Ntpq_rsu = 1;
        return Vtpq_rsu;
    }
    else{
        return Vtpq_rsu;
    }
}

void flexmat6::update_as_tpq_rus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vq*iNt*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vs*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpq_rus(){
    if(Ntpq_rus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vq*iNp*iNt;
        locations.col(1) = vr + vu*iNr + vs*iNu*iNr;
        Vtpq_rus = sp_mat(locations.t(), vValues, iNt*iNp*iNq,iNr*iNu*iNs);
        Ntpq_rus = 1;
        return Vtpq_rus;
    }
    else{
        return Vtpq_rus;
    }
}

void flexmat6::update_as_tpq_sru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vq*iNt*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vu*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpq_sru(){
    if(Ntpq_sru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vq*iNp*iNt;
        locations.col(1) = vs + vr*iNs + vu*iNr*iNs;
        Vtpq_sru = sp_mat(locations.t(), vValues, iNt*iNp*iNq,iNs*iNr*iNu);
        Ntpq_sru = 1;
        return Vtpq_sru;
    }
    else{
        return Vtpq_sru;
    }
}

void flexmat6::update_as_tpq_sur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vq*iNt*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vr*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpq_sur(){
    if(Ntpq_sur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vq*iNp*iNt;
        locations.col(1) = vs + vu*iNs + vr*iNu*iNs;
        Vtpq_sur = sp_mat(locations.t(), vValues, iNt*iNp*iNq,iNs*iNu*iNr);
        Ntpq_sur = 1;
        return Vtpq_sur;
    }
    else{
        return Vtpq_sur;
    }
}

void flexmat6::update_as_tpq_urs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vq*iNt*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vs*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpq_urs(){
    if(Ntpq_urs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vq*iNp*iNt;
        locations.col(1) = vu + vr*iNu + vs*iNr*iNu;
        Vtpq_urs = sp_mat(locations.t(), vValues, iNt*iNp*iNq,iNu*iNr*iNs);
        Ntpq_urs = 1;
        return Vtpq_urs;
    }
    else{
        return Vtpq_urs;
    }
}

void flexmat6::update_as_tpq_usr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vq*iNt*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vr*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpq_usr(){
    if(Ntpq_usr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vq*iNp*iNt;
        locations.col(1) = vu + vs*iNu + vr*iNs*iNu;
        Vtpq_usr = sp_mat(locations.t(), vValues, iNt*iNp*iNq,iNu*iNs*iNr);
        Ntpq_usr = 1;
        return Vtpq_usr;
    }
    else{
        return Vtpq_usr;
    }
}

void flexmat6::update_as_tpr_qsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vr*iNt*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vu*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpr_qsu(){
    if(Ntpr_qsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vr*iNp*iNt;
        locations.col(1) = vq + vs*iNq + vu*iNs*iNq;
        Vtpr_qsu = sp_mat(locations.t(), vValues, iNt*iNp*iNr,iNq*iNs*iNu);
        Ntpr_qsu = 1;
        return Vtpr_qsu;
    }
    else{
        return Vtpr_qsu;
    }
}

void flexmat6::update_as_tpr_qus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vr*iNt*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vs*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpr_qus(){
    if(Ntpr_qus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vr*iNp*iNt;
        locations.col(1) = vq + vu*iNq + vs*iNu*iNq;
        Vtpr_qus = sp_mat(locations.t(), vValues, iNt*iNp*iNr,iNq*iNu*iNs);
        Ntpr_qus = 1;
        return Vtpr_qus;
    }
    else{
        return Vtpr_qus;
    }
}

void flexmat6::update_as_tpr_squ(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vr*iNt*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vu*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpr_squ(){
    if(Ntpr_squ == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vr*iNp*iNt;
        locations.col(1) = vs + vq*iNs + vu*iNq*iNs;
        Vtpr_squ = sp_mat(locations.t(), vValues, iNt*iNp*iNr,iNs*iNq*iNu);
        Ntpr_squ = 1;
        return Vtpr_squ;
    }
    else{
        return Vtpr_squ;
    }
}

void flexmat6::update_as_tpr_suq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vr*iNt*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vq*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpr_suq(){
    if(Ntpr_suq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vr*iNp*iNt;
        locations.col(1) = vs + vu*iNs + vq*iNu*iNs;
        Vtpr_suq = sp_mat(locations.t(), vValues, iNt*iNp*iNr,iNs*iNu*iNq);
        Ntpr_suq = 1;
        return Vtpr_suq;
    }
    else{
        return Vtpr_suq;
    }
}

void flexmat6::update_as_tpr_uqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vr*iNt*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vs*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpr_uqs(){
    if(Ntpr_uqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vr*iNp*iNt;
        locations.col(1) = vu + vq*iNu + vs*iNq*iNu;
        Vtpr_uqs = sp_mat(locations.t(), vValues, iNt*iNp*iNr,iNu*iNq*iNs);
        Ntpr_uqs = 1;
        return Vtpr_uqs;
    }
    else{
        return Vtpr_uqs;
    }
}

void flexmat6::update_as_tpr_usq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vr*iNt*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vq*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpr_usq(){
    if(Ntpr_usq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vr*iNp*iNt;
        locations.col(1) = vu + vs*iNu + vq*iNs*iNu;
        Vtpr_usq = sp_mat(locations.t(), vValues, iNt*iNp*iNr,iNu*iNs*iNq);
        Ntpr_usq = 1;
        return Vtpr_usq;
    }
    else{
        return Vtpr_usq;
    }
}

void flexmat6::update_as_tps_qru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vs*iNt*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vu*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tps_qru(){
    if(Ntps_qru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vs*iNp*iNt;
        locations.col(1) = vq + vr*iNq + vu*iNr*iNq;
        Vtps_qru = sp_mat(locations.t(), vValues, iNt*iNp*iNs,iNq*iNr*iNu);
        Ntps_qru = 1;
        return Vtps_qru;
    }
    else{
        return Vtps_qru;
    }
}

void flexmat6::update_as_tps_qur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vs*iNt*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vr*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tps_qur(){
    if(Ntps_qur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vs*iNp*iNt;
        locations.col(1) = vq + vu*iNq + vr*iNu*iNq;
        Vtps_qur = sp_mat(locations.t(), vValues, iNt*iNp*iNs,iNq*iNu*iNr);
        Ntps_qur = 1;
        return Vtps_qur;
    }
    else{
        return Vtps_qur;
    }
}

void flexmat6::update_as_tps_rqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vs*iNt*iNp)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vu*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tps_rqu(){
    if(Ntps_rqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vs*iNp*iNt;
        locations.col(1) = vr + vq*iNr + vu*iNq*iNr;
        Vtps_rqu = sp_mat(locations.t(), vValues, iNt*iNp*iNs,iNr*iNq*iNu);
        Ntps_rqu = 1;
        return Vtps_rqu;
    }
    else{
        return Vtps_rqu;
    }
}

void flexmat6::update_as_tps_ruq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vs*iNt*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vq*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tps_ruq(){
    if(Ntps_ruq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vs*iNp*iNt;
        locations.col(1) = vr + vu*iNr + vq*iNu*iNr;
        Vtps_ruq = sp_mat(locations.t(), vValues, iNt*iNp*iNs,iNr*iNu*iNq);
        Ntps_ruq = 1;
        return Vtps_ruq;
    }
    else{
        return Vtps_ruq;
    }
}

void flexmat6::update_as_tps_uqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vs*iNt*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vr*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tps_uqr(){
    if(Ntps_uqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vs*iNp*iNt;
        locations.col(1) = vu + vq*iNu + vr*iNq*iNu;
        Vtps_uqr = sp_mat(locations.t(), vValues, iNt*iNp*iNs,iNu*iNq*iNr);
        Ntps_uqr = 1;
        return Vtps_uqr;
    }
    else{
        return Vtps_uqr;
    }
}

void flexmat6::update_as_tps_urq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vs*iNt*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vq*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tps_urq(){
    if(Ntps_urq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vs*iNp*iNt;
        locations.col(1) = vu + vr*iNu + vq*iNr*iNu;
        Vtps_urq = sp_mat(locations.t(), vValues, iNt*iNp*iNs,iNu*iNr*iNq);
        Ntps_urq = 1;
        return Vtps_urq;
    }
    else{
        return Vtps_urq;
    }
}

void flexmat6::update_as_tpu_qrs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vu*iNt*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vs*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpu_qrs(){
    if(Ntpu_qrs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vu*iNp*iNt;
        locations.col(1) = vq + vr*iNq + vs*iNr*iNq;
        Vtpu_qrs = sp_mat(locations.t(), vValues, iNt*iNp*iNu,iNq*iNr*iNs);
        Ntpu_qrs = 1;
        return Vtpu_qrs;
    }
    else{
        return Vtpu_qrs;
    }
}

void flexmat6::update_as_tpu_qsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vu*iNt*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vr*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpu_qsr(){
    if(Ntpu_qsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vu*iNp*iNt;
        locations.col(1) = vq + vs*iNq + vr*iNs*iNq;
        Vtpu_qsr = sp_mat(locations.t(), vValues, iNt*iNp*iNu,iNq*iNs*iNr);
        Ntpu_qsr = 1;
        return Vtpu_qsr;
    }
    else{
        return Vtpu_qsr;
    }
}

void flexmat6::update_as_tpu_rqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vu*iNt*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vs*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpu_rqs(){
    if(Ntpu_rqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vu*iNp*iNt;
        locations.col(1) = vr + vq*iNr + vs*iNq*iNr;
        Vtpu_rqs = sp_mat(locations.t(), vValues, iNt*iNp*iNu,iNr*iNq*iNs);
        Ntpu_rqs = 1;
        return Vtpu_rqs;
    }
    else{
        return Vtpu_rqs;
    }
}

void flexmat6::update_as_tpu_rsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vu*iNt*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vq*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpu_rsq(){
    if(Ntpu_rsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vu*iNp*iNt;
        locations.col(1) = vr + vs*iNr + vq*iNs*iNr;
        Vtpu_rsq = sp_mat(locations.t(), vValues, iNt*iNp*iNu,iNr*iNs*iNq);
        Ntpu_rsq = 1;
        return Vtpu_rsq;
    }
    else{
        return Vtpu_rsq;
    }
}

void flexmat6::update_as_tpu_sqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vu*iNt*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vr*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpu_sqr(){
    if(Ntpu_sqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vu*iNp*iNt;
        locations.col(1) = vs + vq*iNs + vr*iNq*iNs;
        Vtpu_sqr = sp_mat(locations.t(), vValues, iNt*iNp*iNu,iNs*iNq*iNr);
        Ntpu_sqr = 1;
        return Vtpu_sqr;
    }
    else{
        return Vtpu_sqr;
    }
}

void flexmat6::update_as_tpu_srq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNt - vu*iNt*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vq*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tpu_srq(){
    if(Ntpu_srq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vp*iNt + vu*iNp*iNt;
        locations.col(1) = vs + vr*iNs + vq*iNr*iNs;
        Vtpu_srq = sp_mat(locations.t(), vValues, iNt*iNp*iNu,iNs*iNr*iNq);
        Ntpu_srq = 1;
        return Vtpu_srq;
    }
    else{
        return Vtpu_srq;
    }
}

void flexmat6::update_as_tqp_rsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vp*iNt*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vu*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqp_rsu(){
    if(Ntqp_rsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vp*iNq*iNt;
        locations.col(1) = vr + vs*iNr + vu*iNs*iNr;
        Vtqp_rsu = sp_mat(locations.t(), vValues, iNt*iNq*iNp,iNr*iNs*iNu);
        Ntqp_rsu = 1;
        return Vtqp_rsu;
    }
    else{
        return Vtqp_rsu;
    }
}

void flexmat6::update_as_tqp_rus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vp*iNt*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vs*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqp_rus(){
    if(Ntqp_rus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vp*iNq*iNt;
        locations.col(1) = vr + vu*iNr + vs*iNu*iNr;
        Vtqp_rus = sp_mat(locations.t(), vValues, iNt*iNq*iNp,iNr*iNu*iNs);
        Ntqp_rus = 1;
        return Vtqp_rus;
    }
    else{
        return Vtqp_rus;
    }
}

void flexmat6::update_as_tqp_sru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vp*iNt*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vu*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqp_sru(){
    if(Ntqp_sru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vp*iNq*iNt;
        locations.col(1) = vs + vr*iNs + vu*iNr*iNs;
        Vtqp_sru = sp_mat(locations.t(), vValues, iNt*iNq*iNp,iNs*iNr*iNu);
        Ntqp_sru = 1;
        return Vtqp_sru;
    }
    else{
        return Vtqp_sru;
    }
}

void flexmat6::update_as_tqp_sur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vp*iNt*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vr*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqp_sur(){
    if(Ntqp_sur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vp*iNq*iNt;
        locations.col(1) = vs + vu*iNs + vr*iNu*iNs;
        Vtqp_sur = sp_mat(locations.t(), vValues, iNt*iNq*iNp,iNs*iNu*iNr);
        Ntqp_sur = 1;
        return Vtqp_sur;
    }
    else{
        return Vtqp_sur;
    }
}

void flexmat6::update_as_tqp_urs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vp*iNt*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vs*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqp_urs(){
    if(Ntqp_urs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vp*iNq*iNt;
        locations.col(1) = vu + vr*iNu + vs*iNr*iNu;
        Vtqp_urs = sp_mat(locations.t(), vValues, iNt*iNq*iNp,iNu*iNr*iNs);
        Ntqp_urs = 1;
        return Vtqp_urs;
    }
    else{
        return Vtqp_urs;
    }
}

void flexmat6::update_as_tqp_usr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vp*iNt*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vr*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqp_usr(){
    if(Ntqp_usr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vp*iNq*iNt;
        locations.col(1) = vu + vs*iNu + vr*iNs*iNu;
        Vtqp_usr = sp_mat(locations.t(), vValues, iNt*iNq*iNp,iNu*iNs*iNr);
        Ntqp_usr = 1;
        return Vtqp_usr;
    }
    else{
        return Vtqp_usr;
    }
}

void flexmat6::update_as_tqr_psu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vr*iNt*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vu*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqr_psu(){
    if(Ntqr_psu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vr*iNq*iNt;
        locations.col(1) = vp + vs*iNp + vu*iNs*iNp;
        Vtqr_psu = sp_mat(locations.t(), vValues, iNt*iNq*iNr,iNp*iNs*iNu);
        Ntqr_psu = 1;
        return Vtqr_psu;
    }
    else{
        return Vtqr_psu;
    }
}

void flexmat6::update_as_tqr_pus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vr*iNt*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vs*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqr_pus(){
    if(Ntqr_pus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vr*iNq*iNt;
        locations.col(1) = vp + vu*iNp + vs*iNu*iNp;
        Vtqr_pus = sp_mat(locations.t(), vValues, iNt*iNq*iNr,iNp*iNu*iNs);
        Ntqr_pus = 1;
        return Vtqr_pus;
    }
    else{
        return Vtqr_pus;
    }
}

void flexmat6::update_as_tqr_spu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vr*iNt*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vu*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqr_spu(){
    if(Ntqr_spu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vr*iNq*iNt;
        locations.col(1) = vs + vp*iNs + vu*iNp*iNs;
        Vtqr_spu = sp_mat(locations.t(), vValues, iNt*iNq*iNr,iNs*iNp*iNu);
        Ntqr_spu = 1;
        return Vtqr_spu;
    }
    else{
        return Vtqr_spu;
    }
}

void flexmat6::update_as_tqr_sup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vr*iNt*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vp*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqr_sup(){
    if(Ntqr_sup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vr*iNq*iNt;
        locations.col(1) = vs + vu*iNs + vp*iNu*iNs;
        Vtqr_sup = sp_mat(locations.t(), vValues, iNt*iNq*iNr,iNs*iNu*iNp);
        Ntqr_sup = 1;
        return Vtqr_sup;
    }
    else{
        return Vtqr_sup;
    }
}

void flexmat6::update_as_tqr_ups(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vr*iNt*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vs*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqr_ups(){
    if(Ntqr_ups == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vr*iNq*iNt;
        locations.col(1) = vu + vp*iNu + vs*iNp*iNu;
        Vtqr_ups = sp_mat(locations.t(), vValues, iNt*iNq*iNr,iNu*iNp*iNs);
        Ntqr_ups = 1;
        return Vtqr_ups;
    }
    else{
        return Vtqr_ups;
    }
}

void flexmat6::update_as_tqr_usp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vr*iNt*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vp*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqr_usp(){
    if(Ntqr_usp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vr*iNq*iNt;
        locations.col(1) = vu + vs*iNu + vp*iNs*iNu;
        Vtqr_usp = sp_mat(locations.t(), vValues, iNt*iNq*iNr,iNu*iNs*iNp);
        Ntqr_usp = 1;
        return Vtqr_usp;
    }
    else{
        return Vtqr_usp;
    }
}

void flexmat6::update_as_tqs_pru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vs*iNt*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vu*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqs_pru(){
    if(Ntqs_pru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vs*iNq*iNt;
        locations.col(1) = vp + vr*iNp + vu*iNr*iNp;
        Vtqs_pru = sp_mat(locations.t(), vValues, iNt*iNq*iNs,iNp*iNr*iNu);
        Ntqs_pru = 1;
        return Vtqs_pru;
    }
    else{
        return Vtqs_pru;
    }
}

void flexmat6::update_as_tqs_pur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vs*iNt*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vr*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqs_pur(){
    if(Ntqs_pur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vs*iNq*iNt;
        locations.col(1) = vp + vu*iNp + vr*iNu*iNp;
        Vtqs_pur = sp_mat(locations.t(), vValues, iNt*iNq*iNs,iNp*iNu*iNr);
        Ntqs_pur = 1;
        return Vtqs_pur;
    }
    else{
        return Vtqs_pur;
    }
}

void flexmat6::update_as_tqs_rpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vs*iNt*iNq)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vu*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqs_rpu(){
    if(Ntqs_rpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vs*iNq*iNt;
        locations.col(1) = vr + vp*iNr + vu*iNp*iNr;
        Vtqs_rpu = sp_mat(locations.t(), vValues, iNt*iNq*iNs,iNr*iNp*iNu);
        Ntqs_rpu = 1;
        return Vtqs_rpu;
    }
    else{
        return Vtqs_rpu;
    }
}

void flexmat6::update_as_tqs_rup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vs*iNt*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vp*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqs_rup(){
    if(Ntqs_rup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vs*iNq*iNt;
        locations.col(1) = vr + vu*iNr + vp*iNu*iNr;
        Vtqs_rup = sp_mat(locations.t(), vValues, iNt*iNq*iNs,iNr*iNu*iNp);
        Ntqs_rup = 1;
        return Vtqs_rup;
    }
    else{
        return Vtqs_rup;
    }
}

void flexmat6::update_as_tqs_upr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vs*iNt*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vr*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqs_upr(){
    if(Ntqs_upr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vs*iNq*iNt;
        locations.col(1) = vu + vp*iNu + vr*iNp*iNu;
        Vtqs_upr = sp_mat(locations.t(), vValues, iNt*iNq*iNs,iNu*iNp*iNr);
        Ntqs_upr = 1;
        return Vtqs_upr;
    }
    else{
        return Vtqs_upr;
    }
}

void flexmat6::update_as_tqs_urp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vs*iNt*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vp*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqs_urp(){
    if(Ntqs_urp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vs*iNq*iNt;
        locations.col(1) = vu + vr*iNu + vp*iNr*iNu;
        Vtqs_urp = sp_mat(locations.t(), vValues, iNt*iNq*iNs,iNu*iNr*iNp);
        Ntqs_urp = 1;
        return Vtqs_urp;
    }
    else{
        return Vtqs_urp;
    }
}

void flexmat6::update_as_tqu_prs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vu*iNt*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vs*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqu_prs(){
    if(Ntqu_prs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vu*iNq*iNt;
        locations.col(1) = vp + vr*iNp + vs*iNr*iNp;
        Vtqu_prs = sp_mat(locations.t(), vValues, iNt*iNq*iNu,iNp*iNr*iNs);
        Ntqu_prs = 1;
        return Vtqu_prs;
    }
    else{
        return Vtqu_prs;
    }
}

void flexmat6::update_as_tqu_psr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vu*iNt*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vr*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqu_psr(){
    if(Ntqu_psr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vu*iNq*iNt;
        locations.col(1) = vp + vs*iNp + vr*iNs*iNp;
        Vtqu_psr = sp_mat(locations.t(), vValues, iNt*iNq*iNu,iNp*iNs*iNr);
        Ntqu_psr = 1;
        return Vtqu_psr;
    }
    else{
        return Vtqu_psr;
    }
}

void flexmat6::update_as_tqu_rps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vu*iNt*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vs*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqu_rps(){
    if(Ntqu_rps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vu*iNq*iNt;
        locations.col(1) = vr + vp*iNr + vs*iNp*iNr;
        Vtqu_rps = sp_mat(locations.t(), vValues, iNt*iNq*iNu,iNr*iNp*iNs);
        Ntqu_rps = 1;
        return Vtqu_rps;
    }
    else{
        return Vtqu_rps;
    }
}

void flexmat6::update_as_tqu_rsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vu*iNt*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vp*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqu_rsp(){
    if(Ntqu_rsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vu*iNq*iNt;
        locations.col(1) = vr + vs*iNr + vp*iNs*iNr;
        Vtqu_rsp = sp_mat(locations.t(), vValues, iNt*iNq*iNu,iNr*iNs*iNp);
        Ntqu_rsp = 1;
        return Vtqu_rsp;
    }
    else{
        return Vtqu_rsp;
    }
}

void flexmat6::update_as_tqu_spr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vu*iNt*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vr*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqu_spr(){
    if(Ntqu_spr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vu*iNq*iNt;
        locations.col(1) = vs + vp*iNs + vr*iNp*iNs;
        Vtqu_spr = sp_mat(locations.t(), vValues, iNt*iNq*iNu,iNs*iNp*iNr);
        Ntqu_spr = 1;
        return Vtqu_spr;
    }
    else{
        return Vtqu_spr;
    }
}

void flexmat6::update_as_tqu_srp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNt - vu*iNt*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vp*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tqu_srp(){
    if(Ntqu_srp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vq*iNt + vu*iNq*iNt;
        locations.col(1) = vs + vr*iNs + vp*iNr*iNs;
        Vtqu_srp = sp_mat(locations.t(), vValues, iNt*iNq*iNu,iNs*iNr*iNp);
        Ntqu_srp = 1;
        return Vtqu_srp;
    }
    else{
        return Vtqu_srp;
    }
}

void flexmat6::update_as_trp_qsu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vp*iNt*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vu*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trp_qsu(){
    if(Ntrp_qsu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vp*iNr*iNt;
        locations.col(1) = vq + vs*iNq + vu*iNs*iNq;
        Vtrp_qsu = sp_mat(locations.t(), vValues, iNt*iNr*iNp,iNq*iNs*iNu);
        Ntrp_qsu = 1;
        return Vtrp_qsu;
    }
    else{
        return Vtrp_qsu;
    }
}

void flexmat6::update_as_trp_qus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vp*iNt*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vs*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trp_qus(){
    if(Ntrp_qus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vp*iNr*iNt;
        locations.col(1) = vq + vu*iNq + vs*iNu*iNq;
        Vtrp_qus = sp_mat(locations.t(), vValues, iNt*iNr*iNp,iNq*iNu*iNs);
        Ntrp_qus = 1;
        return Vtrp_qus;
    }
    else{
        return Vtrp_qus;
    }
}

void flexmat6::update_as_trp_squ(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vp*iNt*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vu*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trp_squ(){
    if(Ntrp_squ == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vp*iNr*iNt;
        locations.col(1) = vs + vq*iNs + vu*iNq*iNs;
        Vtrp_squ = sp_mat(locations.t(), vValues, iNt*iNr*iNp,iNs*iNq*iNu);
        Ntrp_squ = 1;
        return Vtrp_squ;
    }
    else{
        return Vtrp_squ;
    }
}

void flexmat6::update_as_trp_suq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vp*iNt*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vq*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trp_suq(){
    if(Ntrp_suq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vp*iNr*iNt;
        locations.col(1) = vs + vu*iNs + vq*iNu*iNs;
        Vtrp_suq = sp_mat(locations.t(), vValues, iNt*iNr*iNp,iNs*iNu*iNq);
        Ntrp_suq = 1;
        return Vtrp_suq;
    }
    else{
        return Vtrp_suq;
    }
}

void flexmat6::update_as_trp_uqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vp*iNt*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vs*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trp_uqs(){
    if(Ntrp_uqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vp*iNr*iNt;
        locations.col(1) = vu + vq*iNu + vs*iNq*iNu;
        Vtrp_uqs = sp_mat(locations.t(), vValues, iNt*iNr*iNp,iNu*iNq*iNs);
        Ntrp_uqs = 1;
        return Vtrp_uqs;
    }
    else{
        return Vtrp_uqs;
    }
}

void flexmat6::update_as_trp_usq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vp*iNt*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vq*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trp_usq(){
    if(Ntrp_usq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vp*iNr*iNt;
        locations.col(1) = vu + vs*iNu + vq*iNs*iNu;
        Vtrp_usq = sp_mat(locations.t(), vValues, iNt*iNr*iNp,iNu*iNs*iNq);
        Ntrp_usq = 1;
        return Vtrp_usq;
    }
    else{
        return Vtrp_usq;
    }
}

void flexmat6::update_as_trq_psu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vq*iNt*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vu*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trq_psu(){
    if(Ntrq_psu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vq*iNr*iNt;
        locations.col(1) = vp + vs*iNp + vu*iNs*iNp;
        Vtrq_psu = sp_mat(locations.t(), vValues, iNt*iNr*iNq,iNp*iNs*iNu);
        Ntrq_psu = 1;
        return Vtrq_psu;
    }
    else{
        return Vtrq_psu;
    }
}

void flexmat6::update_as_trq_pus(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vq*iNt*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vs*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trq_pus(){
    if(Ntrq_pus == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vq*iNr*iNt;
        locations.col(1) = vp + vu*iNp + vs*iNu*iNp;
        Vtrq_pus = sp_mat(locations.t(), vValues, iNt*iNr*iNq,iNp*iNu*iNs);
        Ntrq_pus = 1;
        return Vtrq_pus;
    }
    else{
        return Vtrq_pus;
    }
}

void flexmat6::update_as_trq_spu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vq*iNt*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vu*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trq_spu(){
    if(Ntrq_spu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vq*iNr*iNt;
        locations.col(1) = vs + vp*iNs + vu*iNp*iNs;
        Vtrq_spu = sp_mat(locations.t(), vValues, iNt*iNr*iNq,iNs*iNp*iNu);
        Ntrq_spu = 1;
        return Vtrq_spu;
    }
    else{
        return Vtrq_spu;
    }
}

void flexmat6::update_as_trq_sup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vq*iNt*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNu)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vu*iNs - vp*iNs*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trq_sup(){
    if(Ntrq_sup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vq*iNr*iNt;
        locations.col(1) = vs + vu*iNs + vp*iNu*iNs;
        Vtrq_sup = sp_mat(locations.t(), vValues, iNt*iNr*iNq,iNs*iNu*iNp);
        Ntrq_sup = 1;
        return Vtrq_sup;
    }
    else{
        return Vtrq_sup;
    }
}

void flexmat6::update_as_trq_ups(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vq*iNt*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vs*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trq_ups(){
    if(Ntrq_ups == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vq*iNr*iNt;
        locations.col(1) = vu + vp*iNu + vs*iNp*iNu;
        Vtrq_ups = sp_mat(locations.t(), vValues, iNt*iNr*iNq,iNu*iNp*iNs);
        Ntrq_ups = 1;
        return Vtrq_ups;
    }
    else{
        return Vtrq_ups;
    }
}

void flexmat6::update_as_trq_usp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vq*iNt*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vs*iNu - vp*iNu*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trq_usp(){
    if(Ntrq_usp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vq*iNr*iNt;
        locations.col(1) = vu + vs*iNu + vp*iNs*iNu;
        Vtrq_usp = sp_mat(locations.t(), vValues, iNt*iNr*iNq,iNu*iNs*iNp);
        Ntrq_usp = 1;
        return Vtrq_usp;
    }
    else{
        return Vtrq_usp;
    }
}

void flexmat6::update_as_trs_pqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vs*iNt*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vu*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trs_pqu(){
    if(Ntrs_pqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vs*iNr*iNt;
        locations.col(1) = vp + vq*iNp + vu*iNq*iNp;
        Vtrs_pqu = sp_mat(locations.t(), vValues, iNt*iNr*iNs,iNp*iNq*iNu);
        Ntrs_pqu = 1;
        return Vtrs_pqu;
    }
    else{
        return Vtrs_pqu;
    }
}

void flexmat6::update_as_trs_puq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vs*iNt*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vq*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trs_puq(){
    if(Ntrs_puq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vs*iNr*iNt;
        locations.col(1) = vp + vu*iNp + vq*iNu*iNp;
        Vtrs_puq = sp_mat(locations.t(), vValues, iNt*iNr*iNs,iNp*iNu*iNq);
        Ntrs_puq = 1;
        return Vtrs_puq;
    }
    else{
        return Vtrs_puq;
    }
}

void flexmat6::update_as_trs_qpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vs*iNt*iNr)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vu*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trs_qpu(){
    if(Ntrs_qpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vs*iNr*iNt;
        locations.col(1) = vq + vp*iNq + vu*iNp*iNq;
        Vtrs_qpu = sp_mat(locations.t(), vValues, iNt*iNr*iNs,iNq*iNp*iNu);
        Ntrs_qpu = 1;
        return Vtrs_qpu;
    }
    else{
        return Vtrs_qpu;
    }
}

void flexmat6::update_as_trs_qup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vs*iNt*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vp*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trs_qup(){
    if(Ntrs_qup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vs*iNr*iNt;
        locations.col(1) = vq + vu*iNq + vp*iNu*iNq;
        Vtrs_qup = sp_mat(locations.t(), vValues, iNt*iNr*iNs,iNq*iNu*iNp);
        Ntrs_qup = 1;
        return Vtrs_qup;
    }
    else{
        return Vtrs_qup;
    }
}

void flexmat6::update_as_trs_upq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vs*iNt*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vq*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trs_upq(){
    if(Ntrs_upq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vs*iNr*iNt;
        locations.col(1) = vu + vp*iNu + vq*iNp*iNu;
        Vtrs_upq = sp_mat(locations.t(), vValues, iNt*iNr*iNs,iNu*iNp*iNq);
        Ntrs_upq = 1;
        return Vtrs_upq;
    }
    else{
        return Vtrs_upq;
    }
}

void flexmat6::update_as_trs_uqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vs*iNt*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vp*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::trs_uqp(){
    if(Ntrs_uqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vs*iNr*iNt;
        locations.col(1) = vu + vq*iNu + vp*iNq*iNu;
        Vtrs_uqp = sp_mat(locations.t(), vValues, iNt*iNr*iNs,iNu*iNq*iNp);
        Ntrs_uqp = 1;
        return Vtrs_uqp;
    }
    else{
        return Vtrs_uqp;
    }
}

void flexmat6::update_as_tru_pqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vu*iNt*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vs*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tru_pqs(){
    if(Ntru_pqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vu*iNr*iNt;
        locations.col(1) = vp + vq*iNp + vs*iNq*iNp;
        Vtru_pqs = sp_mat(locations.t(), vValues, iNt*iNr*iNu,iNp*iNq*iNs);
        Ntru_pqs = 1;
        return Vtru_pqs;
    }
    else{
        return Vtru_pqs;
    }
}

void flexmat6::update_as_tru_psq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vu*iNt*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vq*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tru_psq(){
    if(Ntru_psq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vu*iNr*iNt;
        locations.col(1) = vp + vs*iNp + vq*iNs*iNp;
        Vtru_psq = sp_mat(locations.t(), vValues, iNt*iNr*iNu,iNp*iNs*iNq);
        Ntru_psq = 1;
        return Vtru_psq;
    }
    else{
        return Vtru_psq;
    }
}

void flexmat6::update_as_tru_qps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vu*iNt*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vs*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tru_qps(){
    if(Ntru_qps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vu*iNr*iNt;
        locations.col(1) = vq + vp*iNq + vs*iNp*iNq;
        Vtru_qps = sp_mat(locations.t(), vValues, iNt*iNr*iNu,iNq*iNp*iNs);
        Ntru_qps = 1;
        return Vtru_qps;
    }
    else{
        return Vtru_qps;
    }
}

void flexmat6::update_as_tru_qsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vu*iNt*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vp*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tru_qsp(){
    if(Ntru_qsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vu*iNr*iNt;
        locations.col(1) = vq + vs*iNq + vp*iNs*iNq;
        Vtru_qsp = sp_mat(locations.t(), vValues, iNt*iNr*iNu,iNq*iNs*iNp);
        Ntru_qsp = 1;
        return Vtru_qsp;
    }
    else{
        return Vtru_qsp;
    }
}

void flexmat6::update_as_tru_spq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vu*iNt*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vq*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tru_spq(){
    if(Ntru_spq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vu*iNr*iNt;
        locations.col(1) = vs + vp*iNs + vq*iNp*iNs;
        Vtru_spq = sp_mat(locations.t(), vValues, iNt*iNr*iNu,iNs*iNp*iNq);
        Ntru_spq = 1;
        return Vtru_spq;
    }
    else{
        return Vtru_spq;
    }
}

void flexmat6::update_as_tru_sqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNt - vu*iNt*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vp*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tru_sqp(){
    if(Ntru_sqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vr*iNt + vu*iNr*iNt;
        locations.col(1) = vs + vq*iNs + vp*iNq*iNs;
        Vtru_sqp = sp_mat(locations.t(), vValues, iNt*iNr*iNu,iNs*iNq*iNp);
        Ntru_sqp = 1;
        return Vtru_sqp;
    }
    else{
        return Vtru_sqp;
    }
}

void flexmat6::update_as_tsp_qru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vp*iNt*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vu*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsp_qru(){
    if(Ntsp_qru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vp*iNs*iNt;
        locations.col(1) = vq + vr*iNq + vu*iNr*iNq;
        Vtsp_qru = sp_mat(locations.t(), vValues, iNt*iNs*iNp,iNq*iNr*iNu);
        Ntsp_qru = 1;
        return Vtsp_qru;
    }
    else{
        return Vtsp_qru;
    }
}

void flexmat6::update_as_tsp_qur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vp*iNt*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vr*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsp_qur(){
    if(Ntsp_qur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vp*iNs*iNt;
        locations.col(1) = vq + vu*iNq + vr*iNu*iNq;
        Vtsp_qur = sp_mat(locations.t(), vValues, iNt*iNs*iNp,iNq*iNu*iNr);
        Ntsp_qur = 1;
        return Vtsp_qur;
    }
    else{
        return Vtsp_qur;
    }
}

void flexmat6::update_as_tsp_rqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vp*iNt*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vu*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsp_rqu(){
    if(Ntsp_rqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vp*iNs*iNt;
        locations.col(1) = vr + vq*iNr + vu*iNq*iNr;
        Vtsp_rqu = sp_mat(locations.t(), vValues, iNt*iNs*iNp,iNr*iNq*iNu);
        Ntsp_rqu = 1;
        return Vtsp_rqu;
    }
    else{
        return Vtsp_rqu;
    }
}

void flexmat6::update_as_tsp_ruq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vp*iNt*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vq*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsp_ruq(){
    if(Ntsp_ruq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vp*iNs*iNt;
        locations.col(1) = vr + vu*iNr + vq*iNu*iNr;
        Vtsp_ruq = sp_mat(locations.t(), vValues, iNt*iNs*iNp,iNr*iNu*iNq);
        Ntsp_ruq = 1;
        return Vtsp_ruq;
    }
    else{
        return Vtsp_ruq;
    }
}

void flexmat6::update_as_tsp_uqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vp*iNt*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vr*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsp_uqr(){
    if(Ntsp_uqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vp*iNs*iNt;
        locations.col(1) = vu + vq*iNu + vr*iNq*iNu;
        Vtsp_uqr = sp_mat(locations.t(), vValues, iNt*iNs*iNp,iNu*iNq*iNr);
        Ntsp_uqr = 1;
        return Vtsp_uqr;
    }
    else{
        return Vtsp_uqr;
    }
}

void flexmat6::update_as_tsp_urq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vp*iNt*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vq*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsp_urq(){
    if(Ntsp_urq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vp*iNs*iNt;
        locations.col(1) = vu + vr*iNu + vq*iNr*iNu;
        Vtsp_urq = sp_mat(locations.t(), vValues, iNt*iNs*iNp,iNu*iNr*iNq);
        Ntsp_urq = 1;
        return Vtsp_urq;
    }
    else{
        return Vtsp_urq;
    }
}

void flexmat6::update_as_tsq_pru(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vq*iNt*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vu*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsq_pru(){
    if(Ntsq_pru == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vq*iNs*iNt;
        locations.col(1) = vp + vr*iNp + vu*iNr*iNp;
        Vtsq_pru = sp_mat(locations.t(), vValues, iNt*iNs*iNq,iNp*iNr*iNu);
        Ntsq_pru = 1;
        return Vtsq_pru;
    }
    else{
        return Vtsq_pru;
    }
}

void flexmat6::update_as_tsq_pur(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vq*iNt*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vr*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsq_pur(){
    if(Ntsq_pur == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vq*iNs*iNt;
        locations.col(1) = vp + vu*iNp + vr*iNu*iNp;
        Vtsq_pur = sp_mat(locations.t(), vValues, iNt*iNs*iNq,iNp*iNu*iNr);
        Ntsq_pur = 1;
        return Vtsq_pur;
    }
    else{
        return Vtsq_pur;
    }
}

void flexmat6::update_as_tsq_rpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vq*iNt*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vu*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsq_rpu(){
    if(Ntsq_rpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vq*iNs*iNt;
        locations.col(1) = vr + vp*iNr + vu*iNp*iNr;
        Vtsq_rpu = sp_mat(locations.t(), vValues, iNt*iNs*iNq,iNr*iNp*iNu);
        Ntsq_rpu = 1;
        return Vtsq_rpu;
    }
    else{
        return Vtsq_rpu;
    }
}

void flexmat6::update_as_tsq_rup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vq*iNt*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNu)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vu*iNr - vp*iNr*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsq_rup(){
    if(Ntsq_rup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vq*iNs*iNt;
        locations.col(1) = vr + vu*iNr + vp*iNu*iNr;
        Vtsq_rup = sp_mat(locations.t(), vValues, iNt*iNs*iNq,iNr*iNu*iNp);
        Ntsq_rup = 1;
        return Vtsq_rup;
    }
    else{
        return Vtsq_rup;
    }
}

void flexmat6::update_as_tsq_upr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vq*iNt*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vr*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsq_upr(){
    if(Ntsq_upr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vq*iNs*iNt;
        locations.col(1) = vu + vp*iNu + vr*iNp*iNu;
        Vtsq_upr = sp_mat(locations.t(), vValues, iNt*iNs*iNq,iNu*iNp*iNr);
        Ntsq_upr = 1;
        return Vtsq_upr;
    }
    else{
        return Vtsq_upr;
    }
}

void flexmat6::update_as_tsq_urp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vq*iNt*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vr*iNu - vp*iNu*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsq_urp(){
    if(Ntsq_urp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vq*iNs*iNt;
        locations.col(1) = vu + vr*iNu + vp*iNr*iNu;
        Vtsq_urp = sp_mat(locations.t(), vValues, iNt*iNs*iNq,iNu*iNr*iNp);
        Ntsq_urp = 1;
        return Vtsq_urp;
    }
    else{
        return Vtsq_urp;
    }
}

void flexmat6::update_as_tsr_pqu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vr*iNt*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vu*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsr_pqu(){
    if(Ntsr_pqu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vr*iNs*iNt;
        locations.col(1) = vp + vq*iNp + vu*iNq*iNp;
        Vtsr_pqu = sp_mat(locations.t(), vValues, iNt*iNs*iNr,iNp*iNq*iNu);
        Ntsr_pqu = 1;
        return Vtsr_pqu;
    }
    else{
        return Vtsr_pqu;
    }
}

void flexmat6::update_as_tsr_puq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vr*iNt*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNu)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNp - vq*iNp*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsr_puq(){
    if(Ntsr_puq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vr*iNs*iNt;
        locations.col(1) = vp + vu*iNp + vq*iNu*iNp;
        Vtsr_puq = sp_mat(locations.t(), vValues, iNt*iNs*iNr,iNp*iNu*iNq);
        Ntsr_puq = 1;
        return Vtsr_puq;
    }
    else{
        return Vtsr_puq;
    }
}

void flexmat6::update_as_tsr_qpu(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vr*iNt*iNs)));
    vu = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vu*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vu*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsr_qpu(){
    if(Ntsr_qpu == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vr*iNs*iNt;
        locations.col(1) = vq + vp*iNq + vu*iNp*iNq;
        Vtsr_qpu = sp_mat(locations.t(), vValues, iNt*iNs*iNr,iNq*iNp*iNu);
        Ntsr_qpu = 1;
        return Vtsr_qpu;
    }
    else{
        return Vtsr_qpu;
    }
}

void flexmat6::update_as_tsr_qup(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vr*iNt*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNu)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vu*iNq - vp*iNq*iNu)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsr_qup(){
    if(Ntsr_qup == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vr*iNs*iNt;
        locations.col(1) = vq + vu*iNq + vp*iNu*iNq;
        Vtsr_qup = sp_mat(locations.t(), vValues, iNt*iNs*iNr,iNq*iNu*iNp);
        Ntsr_qup = 1;
        return Vtsr_qup;
    }
    else{
        return Vtsr_qup;
    }
}

void flexmat6::update_as_tsr_upq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vr*iNt*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vp*iNu - vq*iNu*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsr_upq(){
    if(Ntsr_upq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vr*iNs*iNt;
        locations.col(1) = vu + vp*iNu + vq*iNp*iNu;
        Vtsr_upq = sp_mat(locations.t(), vValues, iNt*iNs*iNr,iNu*iNp*iNq);
        Ntsr_upq = 1;
        return Vtsr_upq;
    }
    else{
        return Vtsr_upq;
    }
}

void flexmat6::update_as_tsr_uqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vr*iNt*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT1 - vq*iNu - vp*iNu*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsr_uqp(){
    if(Ntsr_uqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vr*iNs*iNt;
        locations.col(1) = vu + vq*iNu + vp*iNq*iNu;
        Vtsr_uqp = sp_mat(locations.t(), vValues, iNt*iNs*iNr,iNu*iNq*iNp);
        Ntsr_uqp = 1;
        return Vtsr_uqp;
    }
    else{
        return Vtsr_uqp;
    }
}

void flexmat6::update_as_tsu_pqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vu*iNt*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vr*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsu_pqr(){
    if(Ntsu_pqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vu*iNs*iNt;
        locations.col(1) = vp + vq*iNp + vr*iNq*iNp;
        Vtsu_pqr = sp_mat(locations.t(), vValues, iNt*iNs*iNu,iNp*iNq*iNr);
        Ntsu_pqr = 1;
        return Vtsu_pqr;
    }
    else{
        return Vtsu_pqr;
    }
}

void flexmat6::update_as_tsu_prq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vu*iNt*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vq*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsu_prq(){
    if(Ntsu_prq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vu*iNs*iNt;
        locations.col(1) = vp + vr*iNp + vq*iNr*iNp;
        Vtsu_prq = sp_mat(locations.t(), vValues, iNt*iNs*iNu,iNp*iNr*iNq);
        Ntsu_prq = 1;
        return Vtsu_prq;
    }
    else{
        return Vtsu_prq;
    }
}

void flexmat6::update_as_tsu_qpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vu*iNt*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vr*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsu_qpr(){
    if(Ntsu_qpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vu*iNs*iNt;
        locations.col(1) = vq + vp*iNq + vr*iNp*iNq;
        Vtsu_qpr = sp_mat(locations.t(), vValues, iNt*iNs*iNu,iNq*iNp*iNr);
        Ntsu_qpr = 1;
        return Vtsu_qpr;
    }
    else{
        return Vtsu_qpr;
    }
}

void flexmat6::update_as_tsu_qrp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vu*iNt*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vp*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsu_qrp(){
    if(Ntsu_qrp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vu*iNs*iNt;
        locations.col(1) = vq + vr*iNq + vp*iNr*iNq;
        Vtsu_qrp = sp_mat(locations.t(), vValues, iNt*iNs*iNu,iNq*iNr*iNp);
        Ntsu_qrp = 1;
        return Vtsu_qrp;
    }
    else{
        return Vtsu_qrp;
    }
}

void flexmat6::update_as_tsu_rpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vu*iNt*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vq*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsu_rpq(){
    if(Ntsu_rpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vu*iNs*iNt;
        locations.col(1) = vr + vp*iNr + vq*iNp*iNr;
        Vtsu_rpq = sp_mat(locations.t(), vValues, iNt*iNs*iNu,iNr*iNp*iNq);
        Ntsu_rpq = 1;
        return Vtsu_rpq;
    }
    else{
        return Vtsu_rpq;
    }
}

void flexmat6::update_as_tsu_rqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vu = conv_to<uvec>::from(floor( H.vT0/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vu*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNt - vu*iNt*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vp*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tsu_rqp(){
    if(Ntsu_rqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vs*iNt + vu*iNs*iNt;
        locations.col(1) = vr + vq*iNr + vp*iNq*iNr;
        Vtsu_rqp = sp_mat(locations.t(), vValues, iNt*iNs*iNu,iNr*iNq*iNp);
        Ntsu_rqp = 1;
        return Vtsu_rqp;
    }
    else{
        return Vtsu_rqp;
    }
}

void flexmat6::update_as_tup_qrs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vp*iNt*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vs*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tup_qrs(){
    if(Ntup_qrs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vp*iNu*iNt;
        locations.col(1) = vq + vr*iNq + vs*iNr*iNq;
        Vtup_qrs = sp_mat(locations.t(), vValues, iNt*iNu*iNp,iNq*iNr*iNs);
        Ntup_qrs = 1;
        return Vtup_qrs;
    }
    else{
        return Vtup_qrs;
    }
}

void flexmat6::update_as_tup_qsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vp*iNt*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vr*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tup_qsr(){
    if(Ntup_qsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vp*iNu*iNt;
        locations.col(1) = vq + vs*iNq + vr*iNs*iNq;
        Vtup_qsr = sp_mat(locations.t(), vValues, iNt*iNu*iNp,iNq*iNs*iNr);
        Ntup_qsr = 1;
        return Vtup_qsr;
    }
    else{
        return Vtup_qsr;
    }
}

void flexmat6::update_as_tup_rqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vp*iNt*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vs*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tup_rqs(){
    if(Ntup_rqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vp*iNu*iNt;
        locations.col(1) = vr + vq*iNr + vs*iNq*iNr;
        Vtup_rqs = sp_mat(locations.t(), vValues, iNt*iNu*iNp,iNr*iNq*iNs);
        Ntup_rqs = 1;
        return Vtup_rqs;
    }
    else{
        return Vtup_rqs;
    }
}

void flexmat6::update_as_tup_rsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vp*iNt*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vq*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tup_rsq(){
    if(Ntup_rsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vp*iNu*iNt;
        locations.col(1) = vr + vs*iNr + vq*iNs*iNr;
        Vtup_rsq = sp_mat(locations.t(), vValues, iNt*iNu*iNp,iNr*iNs*iNq);
        Ntup_rsq = 1;
        return Vtup_rsq;
    }
    else{
        return Vtup_rsq;
    }
}

void flexmat6::update_as_tup_sqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vp*iNt*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vr*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tup_sqr(){
    if(Ntup_sqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vp*iNu*iNt;
        locations.col(1) = vs + vq*iNs + vr*iNq*iNs;
        Vtup_sqr = sp_mat(locations.t(), vValues, iNt*iNu*iNp,iNs*iNq*iNr);
        Ntup_sqr = 1;
        return Vtup_sqr;
    }
    else{
        return Vtup_sqr;
    }
}

void flexmat6::update_as_tup_srq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vp*iNt*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vq*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tup_srq(){
    if(Ntup_srq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vp*iNu*iNt;
        locations.col(1) = vs + vr*iNs + vq*iNr*iNs;
        Vtup_srq = sp_mat(locations.t(), vValues, iNt*iNu*iNp,iNs*iNr*iNq);
        Ntup_srq = 1;
        return Vtup_srq;
    }
    else{
        return Vtup_srq;
    }
}

void flexmat6::update_as_tuq_prs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vq*iNt*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vs*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tuq_prs(){
    if(Ntuq_prs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vq*iNu*iNt;
        locations.col(1) = vp + vr*iNp + vs*iNr*iNp;
        Vtuq_prs = sp_mat(locations.t(), vValues, iNt*iNu*iNq,iNp*iNr*iNs);
        Ntuq_prs = 1;
        return Vtuq_prs;
    }
    else{
        return Vtuq_prs;
    }
}

void flexmat6::update_as_tuq_psr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vq*iNt*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vr*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tuq_psr(){
    if(Ntuq_psr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vq*iNu*iNt;
        locations.col(1) = vp + vs*iNp + vr*iNs*iNp;
        Vtuq_psr = sp_mat(locations.t(), vValues, iNt*iNu*iNq,iNp*iNs*iNr);
        Ntuq_psr = 1;
        return Vtuq_psr;
    }
    else{
        return Vtuq_psr;
    }
}

void flexmat6::update_as_tuq_rps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vq*iNt*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vs*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tuq_rps(){
    if(Ntuq_rps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vq*iNu*iNt;
        locations.col(1) = vr + vp*iNr + vs*iNp*iNr;
        Vtuq_rps = sp_mat(locations.t(), vValues, iNt*iNu*iNq,iNr*iNp*iNs);
        Ntuq_rps = 1;
        return Vtuq_rps;
    }
    else{
        return Vtuq_rps;
    }
}

void flexmat6::update_as_tuq_rsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vq*iNt*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vp*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tuq_rsp(){
    if(Ntuq_rsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vq*iNu*iNt;
        locations.col(1) = vr + vs*iNr + vp*iNs*iNr;
        Vtuq_rsp = sp_mat(locations.t(), vValues, iNt*iNu*iNq,iNr*iNs*iNp);
        Ntuq_rsp = 1;
        return Vtuq_rsp;
    }
    else{
        return Vtuq_rsp;
    }
}

void flexmat6::update_as_tuq_spr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vq*iNt*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vr*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tuq_spr(){
    if(Ntuq_spr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vq*iNu*iNt;
        locations.col(1) = vs + vp*iNs + vr*iNp*iNs;
        Vtuq_spr = sp_mat(locations.t(), vValues, iNt*iNu*iNq,iNs*iNp*iNr);
        Ntuq_spr = 1;
        return Vtuq_spr;
    }
    else{
        return Vtuq_spr;
    }
}

void flexmat6::update_as_tuq_srp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vq*iNt*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vp*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tuq_srp(){
    if(Ntuq_srp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vq*iNu*iNt;
        locations.col(1) = vs + vr*iNs + vp*iNr*iNs;
        Vtuq_srp = sp_mat(locations.t(), vValues, iNt*iNu*iNq,iNs*iNr*iNp);
        Ntuq_srp = 1;
        return Vtuq_srp;
    }
    else{
        return Vtuq_srp;
    }
}

void flexmat6::update_as_tur_pqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vr*iNt*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vs*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tur_pqs(){
    if(Ntur_pqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vr*iNu*iNt;
        locations.col(1) = vp + vq*iNp + vs*iNq*iNp;
        Vtur_pqs = sp_mat(locations.t(), vValues, iNt*iNu*iNr,iNp*iNq*iNs);
        Ntur_pqs = 1;
        return Vtur_pqs;
    }
    else{
        return Vtur_pqs;
    }
}

void flexmat6::update_as_tur_psq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vr*iNt*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vq*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tur_psq(){
    if(Ntur_psq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vr*iNu*iNt;
        locations.col(1) = vp + vs*iNp + vq*iNs*iNp;
        Vtur_psq = sp_mat(locations.t(), vValues, iNt*iNu*iNr,iNp*iNs*iNq);
        Ntur_psq = 1;
        return Vtur_psq;
    }
    else{
        return Vtur_psq;
    }
}

void flexmat6::update_as_tur_qps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vr*iNt*iNu)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vs*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tur_qps(){
    if(Ntur_qps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vr*iNu*iNt;
        locations.col(1) = vq + vp*iNq + vs*iNp*iNq;
        Vtur_qps = sp_mat(locations.t(), vValues, iNt*iNu*iNr,iNq*iNp*iNs);
        Ntur_qps = 1;
        return Vtur_qps;
    }
    else{
        return Vtur_qps;
    }
}

void flexmat6::update_as_tur_qsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vr*iNt*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vp*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tur_qsp(){
    if(Ntur_qsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vr*iNu*iNt;
        locations.col(1) = vq + vs*iNq + vp*iNs*iNq;
        Vtur_qsp = sp_mat(locations.t(), vValues, iNt*iNu*iNr,iNq*iNs*iNp);
        Ntur_qsp = 1;
        return Vtur_qsp;
    }
    else{
        return Vtur_qsp;
    }
}

void flexmat6::update_as_tur_spq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vr*iNt*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vq*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tur_spq(){
    if(Ntur_spq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vr*iNu*iNt;
        locations.col(1) = vs + vp*iNs + vq*iNp*iNs;
        Vtur_spq = sp_mat(locations.t(), vValues, iNt*iNu*iNr,iNs*iNp*iNq);
        Ntur_spq = 1;
        return Vtur_spq;
    }
    else{
        return Vtur_spq;
    }
}

void flexmat6::update_as_tur_sqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vr*iNt*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vp*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tur_sqp(){
    if(Ntur_sqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vr*iNu*iNt;
        locations.col(1) = vs + vq*iNs + vp*iNq*iNs;
        Vtur_sqp = sp_mat(locations.t(), vValues, iNt*iNu*iNr,iNs*iNq*iNp);
        Ntur_sqp = 1;
        return Vtur_sqp;
    }
    else{
        return Vtur_sqp;
    }
}

void flexmat6::update_as_tus_pqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vs*iNt*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vr*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tus_pqr(){
    if(Ntus_pqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vs*iNu*iNt;
        locations.col(1) = vp + vq*iNp + vr*iNq*iNp;
        Vtus_pqr = sp_mat(locations.t(), vValues, iNt*iNu*iNs,iNp*iNq*iNr);
        Ntus_pqr = 1;
        return Vtus_pqr;
    }
    else{
        return Vtus_pqr;
    }
}

void flexmat6::update_as_tus_prq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vs*iNt*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vq*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tus_prq(){
    if(Ntus_prq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vs*iNu*iNt;
        locations.col(1) = vp + vr*iNp + vq*iNr*iNp;
        Vtus_prq = sp_mat(locations.t(), vValues, iNt*iNu*iNs,iNp*iNr*iNq);
        Ntus_prq = 1;
        return Vtus_prq;
    }
    else{
        return Vtus_prq;
    }
}

void flexmat6::update_as_tus_qpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vs*iNt*iNu)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vr*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tus_qpr(){
    if(Ntus_qpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vs*iNu*iNt;
        locations.col(1) = vq + vp*iNq + vr*iNp*iNq;
        Vtus_qpr = sp_mat(locations.t(), vValues, iNt*iNu*iNs,iNq*iNp*iNr);
        Ntus_qpr = 1;
        return Vtus_qpr;
    }
    else{
        return Vtus_qpr;
    }
}

void flexmat6::update_as_tus_qrp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vs*iNt*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vp*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tus_qrp(){
    if(Ntus_qrp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vs*iNu*iNt;
        locations.col(1) = vq + vr*iNq + vp*iNr*iNq;
        Vtus_qrp = sp_mat(locations.t(), vValues, iNt*iNu*iNs,iNq*iNr*iNp);
        Ntus_qrp = 1;
        return Vtus_qrp;
    }
    else{
        return Vtus_qrp;
    }
}

void flexmat6::update_as_tus_rpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vs*iNt*iNu)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vq*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tus_rpq(){
    if(Ntus_rpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vs*iNu*iNt;
        locations.col(1) = vr + vp*iNr + vq*iNp*iNr;
        Vtus_rpq = sp_mat(locations.t(), vValues, iNt*iNu*iNs,iNr*iNp*iNq);
        Ntus_rpq = 1;
        return Vtus_rpq;
    }
    else{
        return Vtus_rpq;
    }
}

void flexmat6::update_as_tus_rqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNt*iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNt*iNu)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vu*iNt - vs*iNt*iNu)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vp*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::tus_rqp(){
    if(Ntus_rqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vt + vu*iNt + vs*iNu*iNt;
        locations.col(1) = vr + vq*iNr + vp*iNq*iNr;
        Vtus_rqp = sp_mat(locations.t(), vValues, iNt*iNu*iNs,iNr*iNq*iNp);
        Ntus_rqp = 1;
        return Vtus_rqp;
    }
    else{
        return Vtus_rqp;
    }
}

void flexmat6::update_as_upq_rst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vq*iNu*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vt*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upq_rst(){
    if(Nupq_rst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vq*iNp*iNu;
        locations.col(1) = vr + vs*iNr + vt*iNs*iNr;
        Vupq_rst = sp_mat(locations.t(), vValues, iNu*iNp*iNq,iNr*iNs*iNt);
        Nupq_rst = 1;
        return Vupq_rst;
    }
    else{
        return Vupq_rst;
    }
}

void flexmat6::update_as_upq_rts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vq*iNu*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vs*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upq_rts(){
    if(Nupq_rts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vq*iNp*iNu;
        locations.col(1) = vr + vt*iNr + vs*iNt*iNr;
        Vupq_rts = sp_mat(locations.t(), vValues, iNu*iNp*iNq,iNr*iNt*iNs);
        Nupq_rts = 1;
        return Vupq_rts;
    }
    else{
        return Vupq_rts;
    }
}

void flexmat6::update_as_upq_srt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vq*iNu*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vt*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upq_srt(){
    if(Nupq_srt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vq*iNp*iNu;
        locations.col(1) = vs + vr*iNs + vt*iNr*iNs;
        Vupq_srt = sp_mat(locations.t(), vValues, iNu*iNp*iNq,iNs*iNr*iNt);
        Nupq_srt = 1;
        return Vupq_srt;
    }
    else{
        return Vupq_srt;
    }
}

void flexmat6::update_as_upq_str(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vq*iNu*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vr*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upq_str(){
    if(Nupq_str == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vq*iNp*iNu;
        locations.col(1) = vs + vt*iNs + vr*iNt*iNs;
        Vupq_str = sp_mat(locations.t(), vValues, iNu*iNp*iNq,iNs*iNt*iNr);
        Nupq_str = 1;
        return Vupq_str;
    }
    else{
        return Vupq_str;
    }
}

void flexmat6::update_as_upq_trs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vq*iNu*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vs*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upq_trs(){
    if(Nupq_trs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vq*iNp*iNu;
        locations.col(1) = vt + vr*iNt + vs*iNr*iNt;
        Vupq_trs = sp_mat(locations.t(), vValues, iNu*iNp*iNq,iNt*iNr*iNs);
        Nupq_trs = 1;
        return Vupq_trs;
    }
    else{
        return Vupq_trs;
    }
}

void flexmat6::update_as_upq_tsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vq*iNu*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vr*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upq_tsr(){
    if(Nupq_tsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vq*iNp*iNu;
        locations.col(1) = vt + vs*iNt + vr*iNs*iNt;
        Vupq_tsr = sp_mat(locations.t(), vValues, iNu*iNp*iNq,iNt*iNs*iNr);
        Nupq_tsr = 1;
        return Vupq_tsr;
    }
    else{
        return Vupq_tsr;
    }
}

void flexmat6::update_as_upr_qst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vr*iNu*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vt*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upr_qst(){
    if(Nupr_qst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vr*iNp*iNu;
        locations.col(1) = vq + vs*iNq + vt*iNs*iNq;
        Vupr_qst = sp_mat(locations.t(), vValues, iNu*iNp*iNr,iNq*iNs*iNt);
        Nupr_qst = 1;
        return Vupr_qst;
    }
    else{
        return Vupr_qst;
    }
}

void flexmat6::update_as_upr_qts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vr*iNu*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vs*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upr_qts(){
    if(Nupr_qts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vr*iNp*iNu;
        locations.col(1) = vq + vt*iNq + vs*iNt*iNq;
        Vupr_qts = sp_mat(locations.t(), vValues, iNu*iNp*iNr,iNq*iNt*iNs);
        Nupr_qts = 1;
        return Vupr_qts;
    }
    else{
        return Vupr_qts;
    }
}

void flexmat6::update_as_upr_sqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vr*iNu*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vt*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upr_sqt(){
    if(Nupr_sqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vr*iNp*iNu;
        locations.col(1) = vs + vq*iNs + vt*iNq*iNs;
        Vupr_sqt = sp_mat(locations.t(), vValues, iNu*iNp*iNr,iNs*iNq*iNt);
        Nupr_sqt = 1;
        return Vupr_sqt;
    }
    else{
        return Vupr_sqt;
    }
}

void flexmat6::update_as_upr_stq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vr*iNu*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vq*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upr_stq(){
    if(Nupr_stq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vr*iNp*iNu;
        locations.col(1) = vs + vt*iNs + vq*iNt*iNs;
        Vupr_stq = sp_mat(locations.t(), vValues, iNu*iNp*iNr,iNs*iNt*iNq);
        Nupr_stq = 1;
        return Vupr_stq;
    }
    else{
        return Vupr_stq;
    }
}

void flexmat6::update_as_upr_tqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vr*iNu*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vs*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upr_tqs(){
    if(Nupr_tqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vr*iNp*iNu;
        locations.col(1) = vt + vq*iNt + vs*iNq*iNt;
        Vupr_tqs = sp_mat(locations.t(), vValues, iNu*iNp*iNr,iNt*iNq*iNs);
        Nupr_tqs = 1;
        return Vupr_tqs;
    }
    else{
        return Vupr_tqs;
    }
}

void flexmat6::update_as_upr_tsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vr*iNu*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vq*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upr_tsq(){
    if(Nupr_tsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vr*iNp*iNu;
        locations.col(1) = vt + vs*iNt + vq*iNs*iNt;
        Vupr_tsq = sp_mat(locations.t(), vValues, iNu*iNp*iNr,iNt*iNs*iNq);
        Nupr_tsq = 1;
        return Vupr_tsq;
    }
    else{
        return Vupr_tsq;
    }
}

void flexmat6::update_as_ups_qrt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vs*iNu*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vt*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ups_qrt(){
    if(Nups_qrt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vs*iNp*iNu;
        locations.col(1) = vq + vr*iNq + vt*iNr*iNq;
        Vups_qrt = sp_mat(locations.t(), vValues, iNu*iNp*iNs,iNq*iNr*iNt);
        Nups_qrt = 1;
        return Vups_qrt;
    }
    else{
        return Vups_qrt;
    }
}

void flexmat6::update_as_ups_qtr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vs*iNu*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vr*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ups_qtr(){
    if(Nups_qtr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vs*iNp*iNu;
        locations.col(1) = vq + vt*iNq + vr*iNt*iNq;
        Vups_qtr = sp_mat(locations.t(), vValues, iNu*iNp*iNs,iNq*iNt*iNr);
        Nups_qtr = 1;
        return Vups_qtr;
    }
    else{
        return Vups_qtr;
    }
}

void flexmat6::update_as_ups_rqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vs*iNu*iNp)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vt*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ups_rqt(){
    if(Nups_rqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vs*iNp*iNu;
        locations.col(1) = vr + vq*iNr + vt*iNq*iNr;
        Vups_rqt = sp_mat(locations.t(), vValues, iNu*iNp*iNs,iNr*iNq*iNt);
        Nups_rqt = 1;
        return Vups_rqt;
    }
    else{
        return Vups_rqt;
    }
}

void flexmat6::update_as_ups_rtq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vs*iNu*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vq*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ups_rtq(){
    if(Nups_rtq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vs*iNp*iNu;
        locations.col(1) = vr + vt*iNr + vq*iNt*iNr;
        Vups_rtq = sp_mat(locations.t(), vValues, iNu*iNp*iNs,iNr*iNt*iNq);
        Nups_rtq = 1;
        return Vups_rtq;
    }
    else{
        return Vups_rtq;
    }
}

void flexmat6::update_as_ups_tqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vs*iNu*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vr*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ups_tqr(){
    if(Nups_tqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vs*iNp*iNu;
        locations.col(1) = vt + vq*iNt + vr*iNq*iNt;
        Vups_tqr = sp_mat(locations.t(), vValues, iNu*iNp*iNs,iNt*iNq*iNr);
        Nups_tqr = 1;
        return Vups_tqr;
    }
    else{
        return Vups_tqr;
    }
}

void flexmat6::update_as_ups_trq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vs*iNu*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vq*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ups_trq(){
    if(Nups_trq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vs*iNp*iNu;
        locations.col(1) = vt + vr*iNt + vq*iNr*iNt;
        Vups_trq = sp_mat(locations.t(), vValues, iNu*iNp*iNs,iNt*iNr*iNq);
        Nups_trq = 1;
        return Vups_trq;
    }
    else{
        return Vups_trq;
    }
}

void flexmat6::update_as_upt_qrs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vt*iNu*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vs*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upt_qrs(){
    if(Nupt_qrs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vt*iNp*iNu;
        locations.col(1) = vq + vr*iNq + vs*iNr*iNq;
        Vupt_qrs = sp_mat(locations.t(), vValues, iNu*iNp*iNt,iNq*iNr*iNs);
        Nupt_qrs = 1;
        return Vupt_qrs;
    }
    else{
        return Vupt_qrs;
    }
}

void flexmat6::update_as_upt_qsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vt*iNu*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vr*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upt_qsr(){
    if(Nupt_qsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vt*iNp*iNu;
        locations.col(1) = vq + vs*iNq + vr*iNs*iNq;
        Vupt_qsr = sp_mat(locations.t(), vValues, iNu*iNp*iNt,iNq*iNs*iNr);
        Nupt_qsr = 1;
        return Vupt_qsr;
    }
    else{
        return Vupt_qsr;
    }
}

void flexmat6::update_as_upt_rqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vt*iNu*iNp)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vs*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upt_rqs(){
    if(Nupt_rqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vt*iNp*iNu;
        locations.col(1) = vr + vq*iNr + vs*iNq*iNr;
        Vupt_rqs = sp_mat(locations.t(), vValues, iNu*iNp*iNt,iNr*iNq*iNs);
        Nupt_rqs = 1;
        return Vupt_rqs;
    }
    else{
        return Vupt_rqs;
    }
}

void flexmat6::update_as_upt_rsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vt*iNu*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vq*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upt_rsq(){
    if(Nupt_rsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vt*iNp*iNu;
        locations.col(1) = vr + vs*iNr + vq*iNs*iNr;
        Vupt_rsq = sp_mat(locations.t(), vValues, iNu*iNp*iNt,iNr*iNs*iNq);
        Nupt_rsq = 1;
        return Vupt_rsq;
    }
    else{
        return Vupt_rsq;
    }
}

void flexmat6::update_as_upt_sqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vt*iNu*iNp)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vr*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upt_sqr(){
    if(Nupt_sqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vt*iNp*iNu;
        locations.col(1) = vs + vq*iNs + vr*iNq*iNs;
        Vupt_sqr = sp_mat(locations.t(), vValues, iNu*iNp*iNt,iNs*iNq*iNr);
        Nupt_sqr = 1;
        return Vupt_sqr;
    }
    else{
        return Vupt_sqr;
    }
}

void flexmat6::update_as_upt_srq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNp)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vp*iNu - vt*iNu*iNp)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vq*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::upt_srq(){
    if(Nupt_srq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vp*iNu + vt*iNp*iNu;
        locations.col(1) = vs + vr*iNs + vq*iNr*iNs;
        Vupt_srq = sp_mat(locations.t(), vValues, iNu*iNp*iNt,iNs*iNr*iNq);
        Nupt_srq = 1;
        return Vupt_srq;
    }
    else{
        return Vupt_srq;
    }
}

void flexmat6::update_as_uqp_rst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vp*iNu*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vt*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqp_rst(){
    if(Nuqp_rst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vp*iNq*iNu;
        locations.col(1) = vr + vs*iNr + vt*iNs*iNr;
        Vuqp_rst = sp_mat(locations.t(), vValues, iNu*iNq*iNp,iNr*iNs*iNt);
        Nuqp_rst = 1;
        return Vuqp_rst;
    }
    else{
        return Vuqp_rst;
    }
}

void flexmat6::update_as_uqp_rts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vp*iNu*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vs*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqp_rts(){
    if(Nuqp_rts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vp*iNq*iNu;
        locations.col(1) = vr + vt*iNr + vs*iNt*iNr;
        Vuqp_rts = sp_mat(locations.t(), vValues, iNu*iNq*iNp,iNr*iNt*iNs);
        Nuqp_rts = 1;
        return Vuqp_rts;
    }
    else{
        return Vuqp_rts;
    }
}

void flexmat6::update_as_uqp_srt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vp*iNu*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vt*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqp_srt(){
    if(Nuqp_srt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vp*iNq*iNu;
        locations.col(1) = vs + vr*iNs + vt*iNr*iNs;
        Vuqp_srt = sp_mat(locations.t(), vValues, iNu*iNq*iNp,iNs*iNr*iNt);
        Nuqp_srt = 1;
        return Vuqp_srt;
    }
    else{
        return Vuqp_srt;
    }
}

void flexmat6::update_as_uqp_str(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vp*iNu*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vr*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqp_str(){
    if(Nuqp_str == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vp*iNq*iNu;
        locations.col(1) = vs + vt*iNs + vr*iNt*iNs;
        Vuqp_str = sp_mat(locations.t(), vValues, iNu*iNq*iNp,iNs*iNt*iNr);
        Nuqp_str = 1;
        return Vuqp_str;
    }
    else{
        return Vuqp_str;
    }
}

void flexmat6::update_as_uqp_trs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vp*iNu*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vs*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqp_trs(){
    if(Nuqp_trs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vp*iNq*iNu;
        locations.col(1) = vt + vr*iNt + vs*iNr*iNt;
        Vuqp_trs = sp_mat(locations.t(), vValues, iNu*iNq*iNp,iNt*iNr*iNs);
        Nuqp_trs = 1;
        return Vuqp_trs;
    }
    else{
        return Vuqp_trs;
    }
}

void flexmat6::update_as_uqp_tsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vp*iNu*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vr*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqp_tsr(){
    if(Nuqp_tsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vp*iNq*iNu;
        locations.col(1) = vt + vs*iNt + vr*iNs*iNt;
        Vuqp_tsr = sp_mat(locations.t(), vValues, iNu*iNq*iNp,iNt*iNs*iNr);
        Nuqp_tsr = 1;
        return Vuqp_tsr;
    }
    else{
        return Vuqp_tsr;
    }
}

void flexmat6::update_as_uqr_pst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vr*iNu*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vt*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqr_pst(){
    if(Nuqr_pst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vr*iNq*iNu;
        locations.col(1) = vp + vs*iNp + vt*iNs*iNp;
        Vuqr_pst = sp_mat(locations.t(), vValues, iNu*iNq*iNr,iNp*iNs*iNt);
        Nuqr_pst = 1;
        return Vuqr_pst;
    }
    else{
        return Vuqr_pst;
    }
}

void flexmat6::update_as_uqr_pts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vr*iNu*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vs*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqr_pts(){
    if(Nuqr_pts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vr*iNq*iNu;
        locations.col(1) = vp + vt*iNp + vs*iNt*iNp;
        Vuqr_pts = sp_mat(locations.t(), vValues, iNu*iNq*iNr,iNp*iNt*iNs);
        Nuqr_pts = 1;
        return Vuqr_pts;
    }
    else{
        return Vuqr_pts;
    }
}

void flexmat6::update_as_uqr_spt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vr*iNu*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vt*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqr_spt(){
    if(Nuqr_spt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vr*iNq*iNu;
        locations.col(1) = vs + vp*iNs + vt*iNp*iNs;
        Vuqr_spt = sp_mat(locations.t(), vValues, iNu*iNq*iNr,iNs*iNp*iNt);
        Nuqr_spt = 1;
        return Vuqr_spt;
    }
    else{
        return Vuqr_spt;
    }
}

void flexmat6::update_as_uqr_stp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vr*iNu*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vp*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqr_stp(){
    if(Nuqr_stp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vr*iNq*iNu;
        locations.col(1) = vs + vt*iNs + vp*iNt*iNs;
        Vuqr_stp = sp_mat(locations.t(), vValues, iNu*iNq*iNr,iNs*iNt*iNp);
        Nuqr_stp = 1;
        return Vuqr_stp;
    }
    else{
        return Vuqr_stp;
    }
}

void flexmat6::update_as_uqr_tps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vr*iNu*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vs*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqr_tps(){
    if(Nuqr_tps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vr*iNq*iNu;
        locations.col(1) = vt + vp*iNt + vs*iNp*iNt;
        Vuqr_tps = sp_mat(locations.t(), vValues, iNu*iNq*iNr,iNt*iNp*iNs);
        Nuqr_tps = 1;
        return Vuqr_tps;
    }
    else{
        return Vuqr_tps;
    }
}

void flexmat6::update_as_uqr_tsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vr*iNu*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vp*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqr_tsp(){
    if(Nuqr_tsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vr*iNq*iNu;
        locations.col(1) = vt + vs*iNt + vp*iNs*iNt;
        Vuqr_tsp = sp_mat(locations.t(), vValues, iNu*iNq*iNr,iNt*iNs*iNp);
        Nuqr_tsp = 1;
        return Vuqr_tsp;
    }
    else{
        return Vuqr_tsp;
    }
}

void flexmat6::update_as_uqs_prt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vs*iNu*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vt*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqs_prt(){
    if(Nuqs_prt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vs*iNq*iNu;
        locations.col(1) = vp + vr*iNp + vt*iNr*iNp;
        Vuqs_prt = sp_mat(locations.t(), vValues, iNu*iNq*iNs,iNp*iNr*iNt);
        Nuqs_prt = 1;
        return Vuqs_prt;
    }
    else{
        return Vuqs_prt;
    }
}

void flexmat6::update_as_uqs_ptr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vs*iNu*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vr*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqs_ptr(){
    if(Nuqs_ptr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vs*iNq*iNu;
        locations.col(1) = vp + vt*iNp + vr*iNt*iNp;
        Vuqs_ptr = sp_mat(locations.t(), vValues, iNu*iNq*iNs,iNp*iNt*iNr);
        Nuqs_ptr = 1;
        return Vuqs_ptr;
    }
    else{
        return Vuqs_ptr;
    }
}

void flexmat6::update_as_uqs_rpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vs*iNu*iNq)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vt*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqs_rpt(){
    if(Nuqs_rpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vs*iNq*iNu;
        locations.col(1) = vr + vp*iNr + vt*iNp*iNr;
        Vuqs_rpt = sp_mat(locations.t(), vValues, iNu*iNq*iNs,iNr*iNp*iNt);
        Nuqs_rpt = 1;
        return Vuqs_rpt;
    }
    else{
        return Vuqs_rpt;
    }
}

void flexmat6::update_as_uqs_rtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vs*iNu*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vp*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqs_rtp(){
    if(Nuqs_rtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vs*iNq*iNu;
        locations.col(1) = vr + vt*iNr + vp*iNt*iNr;
        Vuqs_rtp = sp_mat(locations.t(), vValues, iNu*iNq*iNs,iNr*iNt*iNp);
        Nuqs_rtp = 1;
        return Vuqs_rtp;
    }
    else{
        return Vuqs_rtp;
    }
}

void flexmat6::update_as_uqs_tpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vs*iNu*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vr*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqs_tpr(){
    if(Nuqs_tpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vs*iNq*iNu;
        locations.col(1) = vt + vp*iNt + vr*iNp*iNt;
        Vuqs_tpr = sp_mat(locations.t(), vValues, iNu*iNq*iNs,iNt*iNp*iNr);
        Nuqs_tpr = 1;
        return Vuqs_tpr;
    }
    else{
        return Vuqs_tpr;
    }
}

void flexmat6::update_as_uqs_trp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vs*iNu*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vp*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqs_trp(){
    if(Nuqs_trp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vs*iNq*iNu;
        locations.col(1) = vt + vr*iNt + vp*iNr*iNt;
        Vuqs_trp = sp_mat(locations.t(), vValues, iNu*iNq*iNs,iNt*iNr*iNp);
        Nuqs_trp = 1;
        return Vuqs_trp;
    }
    else{
        return Vuqs_trp;
    }
}

void flexmat6::update_as_uqt_prs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vt*iNu*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vs*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqt_prs(){
    if(Nuqt_prs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vt*iNq*iNu;
        locations.col(1) = vp + vr*iNp + vs*iNr*iNp;
        Vuqt_prs = sp_mat(locations.t(), vValues, iNu*iNq*iNt,iNp*iNr*iNs);
        Nuqt_prs = 1;
        return Vuqt_prs;
    }
    else{
        return Vuqt_prs;
    }
}

void flexmat6::update_as_uqt_psr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vt*iNu*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vr*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqt_psr(){
    if(Nuqt_psr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vt*iNq*iNu;
        locations.col(1) = vp + vs*iNp + vr*iNs*iNp;
        Vuqt_psr = sp_mat(locations.t(), vValues, iNu*iNq*iNt,iNp*iNs*iNr);
        Nuqt_psr = 1;
        return Vuqt_psr;
    }
    else{
        return Vuqt_psr;
    }
}

void flexmat6::update_as_uqt_rps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vt*iNu*iNq)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vs*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqt_rps(){
    if(Nuqt_rps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vt*iNq*iNu;
        locations.col(1) = vr + vp*iNr + vs*iNp*iNr;
        Vuqt_rps = sp_mat(locations.t(), vValues, iNu*iNq*iNt,iNr*iNp*iNs);
        Nuqt_rps = 1;
        return Vuqt_rps;
    }
    else{
        return Vuqt_rps;
    }
}

void flexmat6::update_as_uqt_rsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vt*iNu*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vp*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqt_rsp(){
    if(Nuqt_rsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vt*iNq*iNu;
        locations.col(1) = vr + vs*iNr + vp*iNs*iNr;
        Vuqt_rsp = sp_mat(locations.t(), vValues, iNu*iNq*iNt,iNr*iNs*iNp);
        Nuqt_rsp = 1;
        return Vuqt_rsp;
    }
    else{
        return Vuqt_rsp;
    }
}

void flexmat6::update_as_uqt_spr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vt*iNu*iNq)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vr*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqt_spr(){
    if(Nuqt_spr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vt*iNq*iNu;
        locations.col(1) = vs + vp*iNs + vr*iNp*iNs;
        Vuqt_spr = sp_mat(locations.t(), vValues, iNu*iNq*iNt,iNs*iNp*iNr);
        Nuqt_spr = 1;
        return Vuqt_spr;
    }
    else{
        return Vuqt_spr;
    }
}

void flexmat6::update_as_uqt_srp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNq)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vq*iNu - vt*iNu*iNq)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vp*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uqt_srp(){
    if(Nuqt_srp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vq*iNu + vt*iNq*iNu;
        locations.col(1) = vs + vr*iNs + vp*iNr*iNs;
        Vuqt_srp = sp_mat(locations.t(), vValues, iNu*iNq*iNt,iNs*iNr*iNp);
        Nuqt_srp = 1;
        return Vuqt_srp;
    }
    else{
        return Vuqt_srp;
    }
}

void flexmat6::update_as_urp_qst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vp*iNu*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vt*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urp_qst(){
    if(Nurp_qst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vp*iNr*iNu;
        locations.col(1) = vq + vs*iNq + vt*iNs*iNq;
        Vurp_qst = sp_mat(locations.t(), vValues, iNu*iNr*iNp,iNq*iNs*iNt);
        Nurp_qst = 1;
        return Vurp_qst;
    }
    else{
        return Vurp_qst;
    }
}

void flexmat6::update_as_urp_qts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vp*iNu*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vs*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urp_qts(){
    if(Nurp_qts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vp*iNr*iNu;
        locations.col(1) = vq + vt*iNq + vs*iNt*iNq;
        Vurp_qts = sp_mat(locations.t(), vValues, iNu*iNr*iNp,iNq*iNt*iNs);
        Nurp_qts = 1;
        return Vurp_qts;
    }
    else{
        return Vurp_qts;
    }
}

void flexmat6::update_as_urp_sqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vp*iNu*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vt*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urp_sqt(){
    if(Nurp_sqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vp*iNr*iNu;
        locations.col(1) = vs + vq*iNs + vt*iNq*iNs;
        Vurp_sqt = sp_mat(locations.t(), vValues, iNu*iNr*iNp,iNs*iNq*iNt);
        Nurp_sqt = 1;
        return Vurp_sqt;
    }
    else{
        return Vurp_sqt;
    }
}

void flexmat6::update_as_urp_stq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vp*iNu*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vq*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urp_stq(){
    if(Nurp_stq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vp*iNr*iNu;
        locations.col(1) = vs + vt*iNs + vq*iNt*iNs;
        Vurp_stq = sp_mat(locations.t(), vValues, iNu*iNr*iNp,iNs*iNt*iNq);
        Nurp_stq = 1;
        return Vurp_stq;
    }
    else{
        return Vurp_stq;
    }
}

void flexmat6::update_as_urp_tqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vp*iNu*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vs*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urp_tqs(){
    if(Nurp_tqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vp*iNr*iNu;
        locations.col(1) = vt + vq*iNt + vs*iNq*iNt;
        Vurp_tqs = sp_mat(locations.t(), vValues, iNu*iNr*iNp,iNt*iNq*iNs);
        Nurp_tqs = 1;
        return Vurp_tqs;
    }
    else{
        return Vurp_tqs;
    }
}

void flexmat6::update_as_urp_tsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vp*iNu*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vq*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urp_tsq(){
    if(Nurp_tsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vp*iNr*iNu;
        locations.col(1) = vt + vs*iNt + vq*iNs*iNt;
        Vurp_tsq = sp_mat(locations.t(), vValues, iNu*iNr*iNp,iNt*iNs*iNq);
        Nurp_tsq = 1;
        return Vurp_tsq;
    }
    else{
        return Vurp_tsq;
    }
}

void flexmat6::update_as_urq_pst(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vq*iNu*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vt*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urq_pst(){
    if(Nurq_pst == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vq*iNr*iNu;
        locations.col(1) = vp + vs*iNp + vt*iNs*iNp;
        Vurq_pst = sp_mat(locations.t(), vValues, iNu*iNr*iNq,iNp*iNs*iNt);
        Nurq_pst = 1;
        return Vurq_pst;
    }
    else{
        return Vurq_pst;
    }
}

void flexmat6::update_as_urq_pts(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vq*iNu*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vs*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urq_pts(){
    if(Nurq_pts == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vq*iNr*iNu;
        locations.col(1) = vp + vt*iNp + vs*iNt*iNp;
        Vurq_pts = sp_mat(locations.t(), vValues, iNu*iNr*iNq,iNp*iNt*iNs);
        Nurq_pts = 1;
        return Vurq_pts;
    }
    else{
        return Vurq_pts;
    }
}

void flexmat6::update_as_urq_spt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vq*iNu*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vt*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urq_spt(){
    if(Nurq_spt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vq*iNr*iNu;
        locations.col(1) = vs + vp*iNs + vt*iNp*iNs;
        Vurq_spt = sp_mat(locations.t(), vValues, iNu*iNr*iNq,iNs*iNp*iNt);
        Nurq_spt = 1;
        return Vurq_spt;
    }
    else{
        return Vurq_spt;
    }
}

void flexmat6::update_as_urq_stp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vq*iNu*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNt)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vt*iNs - vp*iNs*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urq_stp(){
    if(Nurq_stp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vq*iNr*iNu;
        locations.col(1) = vs + vt*iNs + vp*iNt*iNs;
        Vurq_stp = sp_mat(locations.t(), vValues, iNu*iNr*iNq,iNs*iNt*iNp);
        Nurq_stp = 1;
        return Vurq_stp;
    }
    else{
        return Vurq_stp;
    }
}

void flexmat6::update_as_urq_tps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vq*iNu*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vs*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urq_tps(){
    if(Nurq_tps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vq*iNr*iNu;
        locations.col(1) = vt + vp*iNt + vs*iNp*iNt;
        Vurq_tps = sp_mat(locations.t(), vValues, iNu*iNr*iNq,iNt*iNp*iNs);
        Nurq_tps = 1;
        return Vurq_tps;
    }
    else{
        return Vurq_tps;
    }
}

void flexmat6::update_as_urq_tsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vq*iNu*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNs)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vs*iNt - vp*iNt*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urq_tsp(){
    if(Nurq_tsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vq*iNr*iNu;
        locations.col(1) = vt + vs*iNt + vp*iNs*iNt;
        Vurq_tsp = sp_mat(locations.t(), vValues, iNu*iNr*iNq,iNt*iNs*iNp);
        Nurq_tsp = 1;
        return Vurq_tsp;
    }
    else{
        return Vurq_tsp;
    }
}

void flexmat6::update_as_urs_pqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vs*iNu*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vt*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urs_pqt(){
    if(Nurs_pqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vs*iNr*iNu;
        locations.col(1) = vp + vq*iNp + vt*iNq*iNp;
        Vurs_pqt = sp_mat(locations.t(), vValues, iNu*iNr*iNs,iNp*iNq*iNt);
        Nurs_pqt = 1;
        return Vurs_pqt;
    }
    else{
        return Vurs_pqt;
    }
}

void flexmat6::update_as_urs_ptq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vs*iNu*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vq*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urs_ptq(){
    if(Nurs_ptq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vs*iNr*iNu;
        locations.col(1) = vp + vt*iNp + vq*iNt*iNp;
        Vurs_ptq = sp_mat(locations.t(), vValues, iNu*iNr*iNs,iNp*iNt*iNq);
        Nurs_ptq = 1;
        return Vurs_ptq;
    }
    else{
        return Vurs_ptq;
    }
}

void flexmat6::update_as_urs_qpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vs*iNu*iNr)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vt*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urs_qpt(){
    if(Nurs_qpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vs*iNr*iNu;
        locations.col(1) = vq + vp*iNq + vt*iNp*iNq;
        Vurs_qpt = sp_mat(locations.t(), vValues, iNu*iNr*iNs,iNq*iNp*iNt);
        Nurs_qpt = 1;
        return Vurs_qpt;
    }
    else{
        return Vurs_qpt;
    }
}

void flexmat6::update_as_urs_qtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vs*iNu*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vp*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urs_qtp(){
    if(Nurs_qtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vs*iNr*iNu;
        locations.col(1) = vq + vt*iNq + vp*iNt*iNq;
        Vurs_qtp = sp_mat(locations.t(), vValues, iNu*iNr*iNs,iNq*iNt*iNp);
        Nurs_qtp = 1;
        return Vurs_qtp;
    }
    else{
        return Vurs_qtp;
    }
}

void flexmat6::update_as_urs_tpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vs*iNu*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vq*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urs_tpq(){
    if(Nurs_tpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vs*iNr*iNu;
        locations.col(1) = vt + vp*iNt + vq*iNp*iNt;
        Vurs_tpq = sp_mat(locations.t(), vValues, iNu*iNr*iNs,iNt*iNp*iNq);
        Nurs_tpq = 1;
        return Vurs_tpq;
    }
    else{
        return Vurs_tpq;
    }
}

void flexmat6::update_as_urs_tqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vs*iNu*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vp*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urs_tqp(){
    if(Nurs_tqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vs*iNr*iNu;
        locations.col(1) = vt + vq*iNt + vp*iNq*iNt;
        Vurs_tqp = sp_mat(locations.t(), vValues, iNu*iNr*iNs,iNt*iNq*iNp);
        Nurs_tqp = 1;
        return Vurs_tqp;
    }
    else{
        return Vurs_tqp;
    }
}

void flexmat6::update_as_urt_pqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vt*iNu*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vs*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urt_pqs(){
    if(Nurt_pqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vt*iNr*iNu;
        locations.col(1) = vp + vq*iNp + vs*iNq*iNp;
        Vurt_pqs = sp_mat(locations.t(), vValues, iNu*iNr*iNt,iNp*iNq*iNs);
        Nurt_pqs = 1;
        return Vurt_pqs;
    }
    else{
        return Vurt_pqs;
    }
}

void flexmat6::update_as_urt_psq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vt*iNu*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vq*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urt_psq(){
    if(Nurt_psq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vt*iNr*iNu;
        locations.col(1) = vp + vs*iNp + vq*iNs*iNp;
        Vurt_psq = sp_mat(locations.t(), vValues, iNu*iNr*iNt,iNp*iNs*iNq);
        Nurt_psq = 1;
        return Vurt_psq;
    }
    else{
        return Vurt_psq;
    }
}

void flexmat6::update_as_urt_qps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vt*iNu*iNr)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vs*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urt_qps(){
    if(Nurt_qps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vt*iNr*iNu;
        locations.col(1) = vq + vp*iNq + vs*iNp*iNq;
        Vurt_qps = sp_mat(locations.t(), vValues, iNu*iNr*iNt,iNq*iNp*iNs);
        Nurt_qps = 1;
        return Vurt_qps;
    }
    else{
        return Vurt_qps;
    }
}

void flexmat6::update_as_urt_qsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vt*iNu*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vp*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urt_qsp(){
    if(Nurt_qsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vt*iNr*iNu;
        locations.col(1) = vq + vs*iNq + vp*iNs*iNq;
        Vurt_qsp = sp_mat(locations.t(), vValues, iNu*iNr*iNt,iNq*iNs*iNp);
        Nurt_qsp = 1;
        return Vurt_qsp;
    }
    else{
        return Vurt_qsp;
    }
}

void flexmat6::update_as_urt_spq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vt*iNu*iNr)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vq*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urt_spq(){
    if(Nurt_spq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vt*iNr*iNu;
        locations.col(1) = vs + vp*iNs + vq*iNp*iNs;
        Vurt_spq = sp_mat(locations.t(), vValues, iNu*iNr*iNt,iNs*iNp*iNq);
        Nurt_spq = 1;
        return Vurt_spq;
    }
    else{
        return Vurt_spq;
    }
}

void flexmat6::update_as_urt_sqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNr)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vr*iNu - vt*iNu*iNr)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vp*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::urt_sqp(){
    if(Nurt_sqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vr*iNu + vt*iNr*iNu;
        locations.col(1) = vs + vq*iNs + vp*iNq*iNs;
        Vurt_sqp = sp_mat(locations.t(), vValues, iNu*iNr*iNt,iNs*iNq*iNp);
        Nurt_sqp = 1;
        return Vurt_sqp;
    }
    else{
        return Vurt_sqp;
    }
}

void flexmat6::update_as_usp_qrt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vp*iNu*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vt*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usp_qrt(){
    if(Nusp_qrt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vp*iNs*iNu;
        locations.col(1) = vq + vr*iNq + vt*iNr*iNq;
        Vusp_qrt = sp_mat(locations.t(), vValues, iNu*iNs*iNp,iNq*iNr*iNt);
        Nusp_qrt = 1;
        return Vusp_qrt;
    }
    else{
        return Vusp_qrt;
    }
}

void flexmat6::update_as_usp_qtr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vp*iNu*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vr*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usp_qtr(){
    if(Nusp_qtr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vp*iNs*iNu;
        locations.col(1) = vq + vt*iNq + vr*iNt*iNq;
        Vusp_qtr = sp_mat(locations.t(), vValues, iNu*iNs*iNp,iNq*iNt*iNr);
        Nusp_qtr = 1;
        return Vusp_qtr;
    }
    else{
        return Vusp_qtr;
    }
}

void flexmat6::update_as_usp_rqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vp*iNu*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vt*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usp_rqt(){
    if(Nusp_rqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vp*iNs*iNu;
        locations.col(1) = vr + vq*iNr + vt*iNq*iNr;
        Vusp_rqt = sp_mat(locations.t(), vValues, iNu*iNs*iNp,iNr*iNq*iNt);
        Nusp_rqt = 1;
        return Vusp_rqt;
    }
    else{
        return Vusp_rqt;
    }
}

void flexmat6::update_as_usp_rtq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vp*iNu*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vq*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usp_rtq(){
    if(Nusp_rtq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vp*iNs*iNu;
        locations.col(1) = vr + vt*iNr + vq*iNt*iNr;
        Vusp_rtq = sp_mat(locations.t(), vValues, iNu*iNs*iNp,iNr*iNt*iNq);
        Nusp_rtq = 1;
        return Vusp_rtq;
    }
    else{
        return Vusp_rtq;
    }
}

void flexmat6::update_as_usp_tqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vp*iNu*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vr*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usp_tqr(){
    if(Nusp_tqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vp*iNs*iNu;
        locations.col(1) = vt + vq*iNt + vr*iNq*iNt;
        Vusp_tqr = sp_mat(locations.t(), vValues, iNu*iNs*iNp,iNt*iNq*iNr);
        Nusp_tqr = 1;
        return Vusp_tqr;
    }
    else{
        return Vusp_tqr;
    }
}

void flexmat6::update_as_usp_trq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vp*iNu*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vq*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usp_trq(){
    if(Nusp_trq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vp*iNs*iNu;
        locations.col(1) = vt + vr*iNt + vq*iNr*iNt;
        Vusp_trq = sp_mat(locations.t(), vValues, iNu*iNs*iNp,iNt*iNr*iNq);
        Nusp_trq = 1;
        return Vusp_trq;
    }
    else{
        return Vusp_trq;
    }
}

void flexmat6::update_as_usq_prt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vq*iNu*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vt*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usq_prt(){
    if(Nusq_prt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vq*iNs*iNu;
        locations.col(1) = vp + vr*iNp + vt*iNr*iNp;
        Vusq_prt = sp_mat(locations.t(), vValues, iNu*iNs*iNq,iNp*iNr*iNt);
        Nusq_prt = 1;
        return Vusq_prt;
    }
    else{
        return Vusq_prt;
    }
}

void flexmat6::update_as_usq_ptr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vq*iNu*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vr*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usq_ptr(){
    if(Nusq_ptr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vq*iNs*iNu;
        locations.col(1) = vp + vt*iNp + vr*iNt*iNp;
        Vusq_ptr = sp_mat(locations.t(), vValues, iNu*iNs*iNq,iNp*iNt*iNr);
        Nusq_ptr = 1;
        return Vusq_ptr;
    }
    else{
        return Vusq_ptr;
    }
}

void flexmat6::update_as_usq_rpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vq*iNu*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vt*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usq_rpt(){
    if(Nusq_rpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vq*iNs*iNu;
        locations.col(1) = vr + vp*iNr + vt*iNp*iNr;
        Vusq_rpt = sp_mat(locations.t(), vValues, iNu*iNs*iNq,iNr*iNp*iNt);
        Nusq_rpt = 1;
        return Vusq_rpt;
    }
    else{
        return Vusq_rpt;
    }
}

void flexmat6::update_as_usq_rtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vq*iNu*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNt)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vt*iNr - vp*iNr*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usq_rtp(){
    if(Nusq_rtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vq*iNs*iNu;
        locations.col(1) = vr + vt*iNr + vp*iNt*iNr;
        Vusq_rtp = sp_mat(locations.t(), vValues, iNu*iNs*iNq,iNr*iNt*iNp);
        Nusq_rtp = 1;
        return Vusq_rtp;
    }
    else{
        return Vusq_rtp;
    }
}

void flexmat6::update_as_usq_tpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vq*iNu*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vr*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usq_tpr(){
    if(Nusq_tpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vq*iNs*iNu;
        locations.col(1) = vt + vp*iNt + vr*iNp*iNt;
        Vusq_tpr = sp_mat(locations.t(), vValues, iNu*iNs*iNq,iNt*iNp*iNr);
        Nusq_tpr = 1;
        return Vusq_tpr;
    }
    else{
        return Vusq_tpr;
    }
}

void flexmat6::update_as_usq_trp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vq*iNu*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNr)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vr*iNt - vp*iNt*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usq_trp(){
    if(Nusq_trp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vq*iNs*iNu;
        locations.col(1) = vt + vr*iNt + vp*iNr*iNt;
        Vusq_trp = sp_mat(locations.t(), vValues, iNu*iNs*iNq,iNt*iNr*iNp);
        Nusq_trp = 1;
        return Vusq_trp;
    }
    else{
        return Vusq_trp;
    }
}

void flexmat6::update_as_usr_pqt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vr*iNu*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vt*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usr_pqt(){
    if(Nusr_pqt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vr*iNs*iNu;
        locations.col(1) = vp + vq*iNp + vt*iNq*iNp;
        Vusr_pqt = sp_mat(locations.t(), vValues, iNu*iNs*iNr,iNp*iNq*iNt);
        Nusr_pqt = 1;
        return Vusr_pqt;
    }
    else{
        return Vusr_pqt;
    }
}

void flexmat6::update_as_usr_ptq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vr*iNu*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNt)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNp - vq*iNp*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usr_ptq(){
    if(Nusr_ptq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vr*iNs*iNu;
        locations.col(1) = vp + vt*iNp + vq*iNt*iNp;
        Vusr_ptq = sp_mat(locations.t(), vValues, iNu*iNs*iNr,iNp*iNt*iNq);
        Nusr_ptq = 1;
        return Vusr_ptq;
    }
    else{
        return Vusr_ptq;
    }
}

void flexmat6::update_as_usr_qpt(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vr*iNu*iNs)));
    vt = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vt*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vt*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usr_qpt(){
    if(Nusr_qpt == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vr*iNs*iNu;
        locations.col(1) = vq + vp*iNq + vt*iNp*iNq;
        Vusr_qpt = sp_mat(locations.t(), vValues, iNu*iNs*iNr,iNq*iNp*iNt);
        Nusr_qpt = 1;
        return Vusr_qpt;
    }
    else{
        return Vusr_qpt;
    }
}

void flexmat6::update_as_usr_qtp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vr*iNu*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNt)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vt*iNq - vp*iNq*iNt)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usr_qtp(){
    if(Nusr_qtp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vr*iNs*iNu;
        locations.col(1) = vq + vt*iNq + vp*iNt*iNq;
        Vusr_qtp = sp_mat(locations.t(), vValues, iNu*iNs*iNr,iNq*iNt*iNp);
        Nusr_qtp = 1;
        return Vusr_qtp;
    }
    else{
        return Vusr_qtp;
    }
}

void flexmat6::update_as_usr_tpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vr*iNu*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNt*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNt*iNp)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vp*iNt - vq*iNt*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usr_tpq(){
    if(Nusr_tpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vr*iNs*iNu;
        locations.col(1) = vt + vp*iNt + vq*iNp*iNt;
        Vusr_tpq = sp_mat(locations.t(), vValues, iNu*iNs*iNr,iNt*iNp*iNq);
        Nusr_tpq = 1;
        return Vusr_tpq;
    }
    else{
        return Vusr_tpq;
    }
}

void flexmat6::update_as_usr_tqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vr*iNu*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNt*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNt*iNq)/(iNt)));
    vt = conv_to<uvec>::from(floor((H.vT1 - vq*iNt - vp*iNt*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::usr_tqp(){
    if(Nusr_tqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vr*iNs*iNu;
        locations.col(1) = vt + vq*iNt + vp*iNq*iNt;
        Vusr_tqp = sp_mat(locations.t(), vValues, iNu*iNs*iNr,iNt*iNq*iNp);
        Nusr_tqp = 1;
        return Vusr_tqp;
    }
    else{
        return Vusr_tqp;
    }
}

void flexmat6::update_as_ust_pqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vt*iNu*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vr*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ust_pqr(){
    if(Nust_pqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vt*iNs*iNu;
        locations.col(1) = vp + vq*iNp + vr*iNq*iNp;
        Vust_pqr = sp_mat(locations.t(), vValues, iNu*iNs*iNt,iNp*iNq*iNr);
        Nust_pqr = 1;
        return Vust_pqr;
    }
    else{
        return Vust_pqr;
    }
}

void flexmat6::update_as_ust_prq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vt*iNu*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vq*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ust_prq(){
    if(Nust_prq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vt*iNs*iNu;
        locations.col(1) = vp + vr*iNp + vq*iNr*iNp;
        Vust_prq = sp_mat(locations.t(), vValues, iNu*iNs*iNt,iNp*iNr*iNq);
        Nust_prq = 1;
        return Vust_prq;
    }
    else{
        return Vust_prq;
    }
}

void flexmat6::update_as_ust_qpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vt*iNu*iNs)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vr*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ust_qpr(){
    if(Nust_qpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vt*iNs*iNu;
        locations.col(1) = vq + vp*iNq + vr*iNp*iNq;
        Vust_qpr = sp_mat(locations.t(), vValues, iNu*iNs*iNt,iNq*iNp*iNr);
        Nust_qpr = 1;
        return Vust_qpr;
    }
    else{
        return Vust_qpr;
    }
}

void flexmat6::update_as_ust_qrp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vt*iNu*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vp*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ust_qrp(){
    if(Nust_qrp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vt*iNs*iNu;
        locations.col(1) = vq + vr*iNq + vp*iNr*iNq;
        Vust_qrp = sp_mat(locations.t(), vValues, iNu*iNs*iNt,iNq*iNr*iNp);
        Nust_qrp = 1;
        return Vust_qrp;
    }
    else{
        return Vust_qrp;
    }
}

void flexmat6::update_as_ust_rpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vt*iNu*iNs)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vq*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ust_rpq(){
    if(Nust_rpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vt*iNs*iNu;
        locations.col(1) = vr + vp*iNr + vq*iNp*iNr;
        Vust_rpq = sp_mat(locations.t(), vValues, iNu*iNs*iNt,iNr*iNp*iNq);
        Nust_rpq = 1;
        return Vust_rpq;
    }
    else{
        return Vust_rpq;
    }
}

void flexmat6::update_as_ust_rqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vt = conv_to<uvec>::from(floor( H.vT0/(iNu*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT0 - vt*iNu*iNs)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vs*iNu - vt*iNu*iNs)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vp*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::ust_rqp(){
    if(Nust_rqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vs*iNu + vt*iNs*iNu;
        locations.col(1) = vr + vq*iNr + vp*iNq*iNr;
        Vust_rqp = sp_mat(locations.t(), vValues, iNu*iNs*iNt,iNr*iNq*iNp);
        Nust_rqp = 1;
        return Vust_rqp;
    }
    else{
        return Vust_rqp;
    }
}

void flexmat6::update_as_utp_qrs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vp*iNu*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vs*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utp_qrs(){
    if(Nutp_qrs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vp*iNt*iNu;
        locations.col(1) = vq + vr*iNq + vs*iNr*iNq;
        Vutp_qrs = sp_mat(locations.t(), vValues, iNu*iNt*iNp,iNq*iNr*iNs);
        Nutp_qrs = 1;
        return Vutp_qrs;
    }
    else{
        return Vutp_qrs;
    }
}

void flexmat6::update_as_utp_qsr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vp*iNu*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vr*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utp_qsr(){
    if(Nutp_qsr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vp*iNt*iNu;
        locations.col(1) = vq + vs*iNq + vr*iNs*iNq;
        Vutp_qsr = sp_mat(locations.t(), vValues, iNu*iNt*iNp,iNq*iNs*iNr);
        Nutp_qsr = 1;
        return Vutp_qsr;
    }
    else{
        return Vutp_qsr;
    }
}

void flexmat6::update_as_utp_rqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vp*iNu*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vs*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utp_rqs(){
    if(Nutp_rqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vp*iNt*iNu;
        locations.col(1) = vr + vq*iNr + vs*iNq*iNr;
        Vutp_rqs = sp_mat(locations.t(), vValues, iNu*iNt*iNp,iNr*iNq*iNs);
        Nutp_rqs = 1;
        return Vutp_rqs;
    }
    else{
        return Vutp_rqs;
    }
}

void flexmat6::update_as_utp_rsq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vp*iNu*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vq*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utp_rsq(){
    if(Nutp_rsq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vp*iNt*iNu;
        locations.col(1) = vr + vs*iNr + vq*iNs*iNr;
        Vutp_rsq = sp_mat(locations.t(), vValues, iNu*iNt*iNp,iNr*iNs*iNq);
        Nutp_rsq = 1;
        return Vutp_rsq;
    }
    else{
        return Vutp_rsq;
    }
}

void flexmat6::update_as_utp_sqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vp*iNu*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vr*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utp_sqr(){
    if(Nutp_sqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vp*iNt*iNu;
        locations.col(1) = vs + vq*iNs + vr*iNq*iNs;
        Vutp_sqr = sp_mat(locations.t(), vValues, iNu*iNt*iNp,iNs*iNq*iNr);
        Nutp_sqr = 1;
        return Vutp_sqr;
    }
    else{
        return Vutp_sqr;
    }
}

void flexmat6::update_as_utp_srq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vp = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vp*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vp*iNu*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vq*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utp_srq(){
    if(Nutp_srq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vp*iNt*iNu;
        locations.col(1) = vs + vr*iNs + vq*iNr*iNs;
        Vutp_srq = sp_mat(locations.t(), vValues, iNu*iNt*iNp,iNs*iNr*iNq);
        Nutp_srq = 1;
        return Vutp_srq;
    }
    else{
        return Vutp_srq;
    }
}

void flexmat6::update_as_utq_prs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vq*iNu*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vs*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utq_prs(){
    if(Nutq_prs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vq*iNt*iNu;
        locations.col(1) = vp + vr*iNp + vs*iNr*iNp;
        Vutq_prs = sp_mat(locations.t(), vValues, iNu*iNt*iNq,iNp*iNr*iNs);
        Nutq_prs = 1;
        return Vutq_prs;
    }
    else{
        return Vutq_prs;
    }
}

void flexmat6::update_as_utq_psr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vq*iNu*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vr*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utq_psr(){
    if(Nutq_psr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vq*iNt*iNu;
        locations.col(1) = vp + vs*iNp + vr*iNs*iNp;
        Vutq_psr = sp_mat(locations.t(), vValues, iNu*iNt*iNq,iNp*iNs*iNr);
        Nutq_psr = 1;
        return Vutq_psr;
    }
    else{
        return Vutq_psr;
    }
}

void flexmat6::update_as_utq_rps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vq*iNu*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vs*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utq_rps(){
    if(Nutq_rps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vq*iNt*iNu;
        locations.col(1) = vr + vp*iNr + vs*iNp*iNr;
        Vutq_rps = sp_mat(locations.t(), vValues, iNu*iNt*iNq,iNr*iNp*iNs);
        Nutq_rps = 1;
        return Vutq_rps;
    }
    else{
        return Vutq_rps;
    }
}

void flexmat6::update_as_utq_rsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vq*iNu*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNs)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vs*iNr - vp*iNr*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utq_rsp(){
    if(Nutq_rsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vq*iNt*iNu;
        locations.col(1) = vr + vs*iNr + vp*iNs*iNr;
        Vutq_rsp = sp_mat(locations.t(), vValues, iNu*iNt*iNq,iNr*iNs*iNp);
        Nutq_rsp = 1;
        return Vutq_rsp;
    }
    else{
        return Vutq_rsp;
    }
}

void flexmat6::update_as_utq_spr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vq*iNu*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vr*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utq_spr(){
    if(Nutq_spr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vq*iNt*iNu;
        locations.col(1) = vs + vp*iNs + vr*iNp*iNs;
        Vutq_spr = sp_mat(locations.t(), vValues, iNu*iNt*iNq,iNs*iNp*iNr);
        Nutq_spr = 1;
        return Vutq_spr;
    }
    else{
        return Vutq_spr;
    }
}

void flexmat6::update_as_utq_srp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vq = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vq*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vq*iNu*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNr)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vr*iNs - vp*iNs*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utq_srp(){
    if(Nutq_srp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vq*iNt*iNu;
        locations.col(1) = vs + vr*iNs + vp*iNr*iNs;
        Vutq_srp = sp_mat(locations.t(), vValues, iNu*iNt*iNq,iNs*iNr*iNp);
        Nutq_srp = 1;
        return Vutq_srp;
    }
    else{
        return Vutq_srp;
    }
}

void flexmat6::update_as_utr_pqs(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vr*iNu*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vs*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utr_pqs(){
    if(Nutr_pqs == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vr*iNt*iNu;
        locations.col(1) = vp + vq*iNp + vs*iNq*iNp;
        Vutr_pqs = sp_mat(locations.t(), vValues, iNu*iNt*iNr,iNp*iNq*iNs);
        Nutr_pqs = 1;
        return Vutr_pqs;
    }
    else{
        return Vutr_pqs;
    }
}

void flexmat6::update_as_utr_psq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vr*iNu*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNs)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNp - vq*iNp*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utr_psq(){
    if(Nutr_psq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vr*iNt*iNu;
        locations.col(1) = vp + vs*iNp + vq*iNs*iNp;
        Vutr_psq = sp_mat(locations.t(), vValues, iNu*iNt*iNr,iNp*iNs*iNq);
        Nutr_psq = 1;
        return Vutr_psq;
    }
    else{
        return Vutr_psq;
    }
}

void flexmat6::update_as_utr_qps(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vr*iNu*iNt)));
    vs = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vs*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vs*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utr_qps(){
    if(Nutr_qps == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vr*iNt*iNu;
        locations.col(1) = vq + vp*iNq + vs*iNp*iNq;
        Vutr_qps = sp_mat(locations.t(), vValues, iNu*iNt*iNr,iNq*iNp*iNs);
        Nutr_qps = 1;
        return Vutr_qps;
    }
    else{
        return Vutr_qps;
    }
}

void flexmat6::update_as_utr_qsp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vr*iNu*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNs)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vs*iNq - vp*iNq*iNs)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utr_qsp(){
    if(Nutr_qsp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vr*iNt*iNu;
        locations.col(1) = vq + vs*iNq + vp*iNs*iNq;
        Vutr_qsp = sp_mat(locations.t(), vValues, iNu*iNt*iNr,iNq*iNs*iNp);
        Nutr_qsp = 1;
        return Vutr_qsp;
    }
    else{
        return Vutr_qsp;
    }
}

void flexmat6::update_as_utr_spq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vr*iNu*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNs*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNs*iNp)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vp*iNs - vq*iNs*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utr_spq(){
    if(Nutr_spq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vr*iNt*iNu;
        locations.col(1) = vs + vp*iNs + vq*iNp*iNs;
        Vutr_spq = sp_mat(locations.t(), vValues, iNu*iNt*iNr,iNs*iNp*iNq);
        Nutr_spq = 1;
        return Vutr_spq;
    }
    else{
        return Vutr_spq;
    }
}

void flexmat6::update_as_utr_sqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vr = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vr*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vr*iNu*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNs*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNs*iNq)/(iNs)));
    vs = conv_to<uvec>::from(floor((H.vT1 - vq*iNs - vp*iNs*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::utr_sqp(){
    if(Nutr_sqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vr*iNt*iNu;
        locations.col(1) = vs + vq*iNs + vp*iNq*iNs;
        Vutr_sqp = sp_mat(locations.t(), vValues, iNu*iNt*iNr,iNs*iNq*iNp);
        Nutr_sqp = 1;
        return Vutr_sqp;
    }
    else{
        return Vutr_sqp;
    }
}

void flexmat6::update_as_uts_pqr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vs*iNu*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNp*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNp*iNq)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNp - vr*iNp*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uts_pqr(){
    if(Nuts_pqr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vs*iNt*iNu;
        locations.col(1) = vp + vq*iNp + vr*iNq*iNp;
        Vuts_pqr = sp_mat(locations.t(), vValues, iNu*iNt*iNs,iNp*iNq*iNr);
        Nuts_pqr = 1;
        return Vuts_pqr;
    }
    else{
        return Vuts_pqr;
    }
}

void flexmat6::update_as_uts_prq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vs*iNu*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNp*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNp*iNr)/(iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNp - vq*iNp*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uts_prq(){
    if(Nuts_prq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vs*iNt*iNu;
        locations.col(1) = vp + vr*iNp + vq*iNr*iNp;
        Vuts_prq = sp_mat(locations.t(), vValues, iNu*iNt*iNs,iNp*iNr*iNq);
        Nuts_prq = 1;
        return Vuts_prq;
    }
    else{
        return Vuts_prq;
    }
}

void flexmat6::update_as_uts_qpr(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vs*iNu*iNt)));
    vr = conv_to<uvec>::from(floor( H.vT1/(iNq*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vr*iNq*iNp)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNq - vr*iNq*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uts_qpr(){
    if(Nuts_qpr == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vs*iNt*iNu;
        locations.col(1) = vq + vp*iNq + vr*iNp*iNq;
        Vuts_qpr = sp_mat(locations.t(), vValues, iNu*iNt*iNs,iNq*iNp*iNr);
        Nuts_qpr = 1;
        return Vuts_qpr;
    }
    else{
        return Vuts_qpr;
    }
}

void flexmat6::update_as_uts_qrp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vs*iNu*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNq*iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNq*iNr)/(iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vr*iNq - vp*iNq*iNr)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uts_qrp(){
    if(Nuts_qrp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vs*iNt*iNu;
        locations.col(1) = vq + vr*iNq + vp*iNr*iNq;
        Vuts_qrp = sp_mat(locations.t(), vValues, iNu*iNt*iNs,iNq*iNr*iNp);
        Nuts_qrp = 1;
        return Vuts_qrp;
    }
    else{
        return Vuts_qrp;
    }
}

void flexmat6::update_as_uts_rpq(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vs*iNu*iNt)));
    vq = conv_to<uvec>::from(floor( H.vT1/(iNr*iNp)));
    vp = conv_to<uvec>::from(floor((H.vT1 - vq*iNr*iNp)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vp*iNr - vq*iNr*iNp)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uts_rpq(){
    if(Nuts_rpq == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vs*iNt*iNu;
        locations.col(1) = vr + vp*iNr + vq*iNp*iNr;
        Vuts_rpq = sp_mat(locations.t(), vValues, iNu*iNt*iNs,iNr*iNp*iNq);
        Nuts_rpq = 1;
        return Vuts_rpq;
    }
    else{
        return Vuts_rpq;
    }
}

void flexmat6::update_as_uts_rqp(sp_mat spC, int Np, int Nq, int Nr, int Ns, int Nt, int Nu){
    iNp = Np;
    iNq = Nq;
    iNr = Nr;
    iNs = Ns;
    iNt = Nt;
    iNu = Nu;
    unpack_sp_mat H(spC);
    vs = conv_to<uvec>::from(floor( H.vT0/(iNu*iNt)));
    vt = conv_to<uvec>::from(floor((H.vT0 - vs*iNu*iNt)/(iNu)));
    vu = conv_to<uvec>::from(floor((H.vT0 - vt*iNu - vs*iNu*iNt)));
    vp = conv_to<uvec>::from(floor( H.vT1/(iNr*iNq)));
    vq = conv_to<uvec>::from(floor((H.vT1 - vp*iNr*iNq)/(iNr)));
    vr = conv_to<uvec>::from(floor((H.vT1 - vq*iNr - vp*iNr*iNq)));
    vValues = H.vVals; 
    deinit();
}

sp_mat flexmat6::uts_rqp(){
    if(Nuts_rqp == 0){
        locations.set_size(vp.size(), 2);
        locations.col(0) = vu + vt*iNu + vs*iNt*iNu;
        locations.col(1) = vr + vq*iNr + vp*iNq*iNr;
        Vuts_rqp = sp_mat(locations.t(), vValues, iNu*iNt*iNs,iNr*iNq*iNp);
        Nuts_rqp = 1;
        return Vuts_rqp;
    }
    else{
        return Vuts_rqp;
    }
}

