'''
Lep1fatJet2_GenPart_eta
Lep1fatJet2_GenPart_genPartIdxMother
Lep1fatJet2_GenPart_mass
Lep1fatJet2_GenPart_pdgId
Lep1fatJet2_GenPart_phi
Lep1fatJet2_GenPart_pt
Lep1fatJet2_GenPart_status
Lep1fatJet2_GenPart_statusFlags
Lep1fatJet2_nGenPart
'''


matching = '''
int matching[2] = {-99};

float ptgenwl[5][15];
for(size_t ik=0; ik<5;ik++){
    for(size_t ij=0; ij<15;ij++){
        ptgenwl[ik][ij] = -99;
    }
}
int Windex = 0;

float ptgentl[5][13] = {-99.};
for(size_t ik=0; ik<5;ik++){
    for(size_t ij=0; ij<13;ij++){
        ptgentl[ik][ij] = -99;
    }
}
int Tindex = 0;

float ptgengqf[14][4] = {-99.};
for(size_t ik=0; ik<14;ik++){
    for(size_t ij=0; ij<4;ij++){
        ptgengqf[ik][ij] = -99;
    }
}
int gqindex = 0;



std::vector<int> FirstCopy_g;
int PDGIDqg[6] = {1,2,3,4,5,21};

for(size_t ik=0; ik<Lep1fatJet2_GenPart_pt.size();ik++){
    // W
    if ( abs(Lep1fatJet2_GenPart_pdgId[ik]) == 24 ){
        if (not (Lep1fatJet2_GenPart_statusFlags[ik]&(1<<13))) continue; // isLastCopy
        if( Windex < 5 ){

            ptgenwl[Windex][0] = Lep1fatJet2_GenPart_pt[ik];
            ptgenwl[Windex][1] = Lep1fatJet2_GenPart_eta[ik];
            ptgenwl[Windex][2] = Lep1fatJet2_GenPart_phi[ik];
            ptgenwl[Windex][3] = Lep1fatJet2_GenPart_mass[ik];
            ptgenwl[Windex][14] = ik;

            vector<int> W_daughter_index;
            for (size_t id=0; id<Lep1fatJet2_GenPart_pt.size();id++){
                if (Lep1fatJet2_GenPart_genPartIdxMother[id] == ik){
                    W_daughter_index.push_back(id);
                }
            }
            int NW_daughter = W_daughter_index.size();
            if ( NW_daughter == 2){
                ptgenwl[Windex][4] = Lep1fatJet2_GenPart_pt[W_daughter_index[0]];
                ptgenwl[Windex][5] = Lep1fatJet2_GenPart_eta[W_daughter_index[0]];
                ptgenwl[Windex][6] = Lep1fatJet2_GenPart_phi[W_daughter_index[0]];
                ptgenwl[Windex][7] = Lep1fatJet2_GenPart_mass[W_daughter_index[0]];
                ptgenwl[Windex][8] = Lep1fatJet2_GenPart_pdgId[W_daughter_index[0]];

                ptgenwl[Windex][9]  = Lep1fatJet2_GenPart_pt[W_daughter_index[1]];
                ptgenwl[Windex][10] = Lep1fatJet2_GenPart_eta[W_daughter_index[1]];
                ptgenwl[Windex][11] = Lep1fatJet2_GenPart_phi[W_daughter_index[1]];
                ptgenwl[Windex][12] = Lep1fatJet2_GenPart_mass[W_daughter_index[1]];
                ptgenwl[Windex][13] = Lep1fatJet2_GenPart_pdgId[W_daughter_index[1]];
            }
            Windex++;
        }
    }

    // q,g
    for(size_t igq=0; igq<6;igq++){
        int PGDID = PDGIDqg[igq];
        if (abs(Lep1fatJet2_GenPart_pdgId[ik]) == PGDID ){
            if( Lep1fatJet2_GenPart_pt[ik] > 50 ){
                if( gqindex < 15 ){
                    bool overlap = false;
                    int FirstCopy = -99;
                    int LoopID = ik;
                    bool From_WZTop = false;
                    while(abs(Lep1fatJet2_GenPart_pdgId[LoopID]) == PGDID){
                        FirstCopy = LoopID;
                        LoopID = Lep1fatJet2_GenPart_genPartIdxMother[LoopID];
                    }
                    while( LoopID >= 0 ){
                        if( abs(Lep1fatJet2_GenPart_pdgId[LoopID]) == 24 || abs(Lep1fatJet2_GenPart_pdgId[LoopID]) == 23 || abs(Lep1fatJet2_GenPart_pdgId[LoopID]) == 6 ){
                            From_WZTop = true;
                        }
                        LoopID = Lep1fatJet2_GenPart_genPartIdxMother[LoopID];
                    }
                    if(From_WZTop) continue;
                    for(unsigned int inum = 0; inum < FirstCopy_g.size() ; ++inum){
                        if( FirstCopy == FirstCopy_g[inum] ){
                            overlap = true;
                        }
                    }
                    if(overlap) continue;
                    FirstCopy_g.push_back(FirstCopy);
                    ptgengqf[gqindex][0] = Lep1fatJet2_GenPart_pt[FirstCopy];
                    ptgengqf[gqindex][1] = Lep1fatJet2_GenPart_eta[FirstCopy];
                    ptgengqf[gqindex][2] = Lep1fatJet2_GenPart_phi[FirstCopy];
                    ptgengqf[gqindex][3] = Lep1fatJet2_GenPart_mass[FirstCopy];
                    gqindex++;
                }
            }
        }
    }
}

for(size_t ik=0; ik<Lep1fatJet2_GenPart_pt.size();ik++){
    // Top
    if ( abs(Lep1fatJet2_GenPart_pdgId[ik]) == 6 ){
        if( Tindex < 4 ){
            bool lastcopy = true;
            vector<int> T_daughter_index;
            for (size_t id=0; id<Lep1fatJet2_GenPart_pt.size();id++){
                if (Lep1fatJet2_GenPart_genPartIdxMother[id] == ik){
                    T_daughter_index.push_back(id);
                }
            }
            int NT_daughter = T_daughter_index.size();
            for (size_t id=0; id<NT_daughter;id++){
                if(abs(Lep1fatJet2_GenPart_pdgId[T_daughter_index[id]])==6){
                    lastcopy = false;
                }
            }
            if(!lastcopy) continue;

            ptgentl[Tindex][0] = Lep1fatJet2_GenPart_pt[ik];
            ptgentl[Tindex][1] = Lep1fatJet2_GenPart_eta[ik];
            ptgentl[Tindex][2] = Lep1fatJet2_GenPart_phi[ik];
            ptgentl[Tindex][3] = Lep1fatJet2_GenPart_mass[ik];
            
            if ( NT_daughter == 2){
                int W_daughter = -99;
                int bdaughter = -99;
                if( abs(Lep1fatJet2_GenPart_pdgId[T_daughter_index[0]]) == 24 ){
                    W_daughter = T_daughter_index[0];
                }
                if( abs(Lep1fatJet2_GenPart_pdgId[T_daughter_index[1]]) == 24 ){
                    W_daughter = T_daughter_index[1];
                }
                if( abs(Lep1fatJet2_GenPart_pdgId[T_daughter_index[0]]) == 5 ){
                    bdaughter = T_daughter_index[0];
                }
                if( abs(Lep1fatJet2_GenPart_pdgId[T_daughter_index[1]]) == 5 ){
                    bdaughter = T_daughter_index[1];
                }

                for(size_t iw=0; iw<5;iw++){
                    if(ptgenwl[iw][0] > 0){
                        int LoopIDW = ptgenwl[iw][14];
                        bool W_From_T = false;
                        while( LoopIDW > 0 ){
                            if( LoopIDW == W_daughter ){
                                W_From_T = true;
                            }
                            LoopIDW = Lep1fatJet2_GenPart_genPartIdxMother[LoopIDW];
                        }
                        if(W_From_T){
                            ptgentl[Tindex][4] = iw;
                        }
                    }
                }
                
                ptgentl[Tindex][5] = Lep1fatJet2_GenPart_pt[bdaughter];
                ptgentl[Tindex][6] = Lep1fatJet2_GenPart_eta[bdaughter];
                ptgentl[Tindex][7] = Lep1fatJet2_GenPart_phi[bdaughter];
                ptgentl[Tindex][8] = Lep1fatJet2_GenPart_mass[bdaughter];

                ptgentl[Tindex][9] = Lep1fatJet2_GenPart_pt[W_daughter];
                ptgentl[Tindex][10] = Lep1fatJet2_GenPart_eta[W_daughter];
                ptgentl[Tindex][11] = Lep1fatJet2_GenPart_phi[W_daughter];
                ptgentl[Tindex][12] = Lep1fatJet2_GenPart_mass[W_daughter];
            }
            Tindex++;
        }
    }
}

TLorentzVector Genpart, Genpartd1, Genpartd2, Genpartd1d1, Genpartd1d2, fJ1, fJ2, lepton;
 
fJ1.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt, Lep1fatJet2_FatJet_eta, Lep1fatJet2_FatJet_phi, Lep1fatJet2_FatJet_msoftdrop );
fJ2.SetPtEtaPhiM( Lep1fatJet2_FatJet_pt_2, Lep1fatJet2_FatJet_eta_2, Lep1fatJet2_FatJet_phi_2, Lep1fatJet2_FatJet_msoftdrop_2 );
lepton.SetPtEtaPhiE( Lep1fatJet2_LeptonPt, Lep1fatJet2_LeptonEta, Lep1fatJet2_LeptonPhi, Lep1fatJet2_LeptonE );

for(size_t ik=0; ik<5;ik++){
    if(ptgenwl[ik][0] > 0){
        Genpart.SetPtEtaPhiM( ptgenwl[ik][0], ptgenwl[ik][1], ptgenwl[ik][2], ptgenwl[ik][3] );
        Genpartd1.SetPtEtaPhiM( ptgenwl[ik][4], ptgenwl[ik][5], ptgenwl[ik][6], ptgenwl[ik][7] );
        Genpartd2.SetPtEtaPhiM( ptgenwl[ik][9], ptgenwl[ik][10], ptgenwl[ik][11], ptgenwl[ik][12] );
        if( fabs(Genpart.DeltaR(fJ1)) < 0.6 ){
            matching[0] = 1;
            if( fabs(Genpartd1.DeltaR(fJ1)) < 0.8 && fabs(Genpartd2.DeltaR(fJ1)) < 0.8 ){
                matching[0] = 2;
            }
        }
        if( ((fabs(Genpartd1.DeltaR(fJ1)) < 0.8) + (fabs(Genpartd2.DeltaR(fJ1)) < 0.8)) == 1 ){
            matching[0] = 4;
        }
        if( Genpart.DeltaR(fJ2) < 0.6 ){
            matching[1] = 1;
            if( fabs(Genpartd1.DeltaR(fJ2)) < 0.8 && fabs(Genpartd2.DeltaR(fJ2)) < 0.8 ){
                matching[1] = 2;
            }
            
        }
        if( ((fabs(Genpartd1.DeltaR(fJ2)) < 0.8) + (fabs(Genpartd2.DeltaR(fJ2)) < 0.8)) == 1 ){
            matching[1] = 4;
        }
    }
}

for(size_t ik=0; ik<15;ik++){
    if(ptgengqf[ik][0] > 50){
        Genpart.SetPtEtaPhiM( ptgengqf[ik][0], ptgengqf[ik][1], ptgengqf[ik][2], ptgengqf[ik][3] );
        if( fabs(Genpart.DeltaR(fJ1)) < 0.6 ){
            matching[0] = 3;
        }
        if( fabs(Genpart.DeltaR(fJ2)) < 0.6 ){
            matching[1] = 3;
        }
    }
}

for(size_t ik=0; ik<4;ik++){
    if(ptgentl[ik][0] > 0){
        Genpart.SetPtEtaPhiM( ptgentl[ik][0], ptgentl[ik][1], ptgentl[ik][2], ptgentl[ik][3] );
        int W_index = ptgentl[ik][4];
        Genpartd1.SetPtEtaPhiM( ptgenwl[W_index][0], ptgenwl[W_index][1], ptgenwl[W_index][2], ptgenwl[W_index][3] );
        Genpartd2.SetPtEtaPhiM( ptgentl[ik][5], ptgenwl[ik][6], ptgenwl[ik][7], ptgenwl[ik][8] );
        Genpartd1d1.SetPtEtaPhiM( ptgenwl[W_index][4], ptgenwl[W_index][5], ptgenwl[W_index][6], ptgenwl[W_index][7] );
        Genpartd1d2.SetPtEtaPhiM( ptgenwl[W_index][9], ptgenwl[W_index][10], ptgenwl[W_index][11], ptgenwl[W_index][12] );

        if( fabs(Genpart.DeltaR(fJ1)) < 0.6 ){
            if( fabs(Genpartd1d1.DeltaR(fJ1)) < 0.8 && fabs(Genpartd1d2.DeltaR(fJ1)) < 0.8 && fabs(Genpartd2.DeltaR(fJ1)) < 0.8 ){
                matching[0] = 5;
            }
            if( ( (fabs(Genpartd1d1.DeltaR(fJ1)) < 0.8)+(fabs(Genpartd1d2.DeltaR(fJ1)) < 0.8) )==1 && fabs(Genpartd2.DeltaR(fJ1)) < 0.8 ){
                matching[0] = 6;
            }
        }

        if( fabs(Genpart.DeltaR(fJ2)) < 0.6 ){
            if( fabs(Genpartd1d1.DeltaR(fJ2)) < 0.8 && fabs(Genpartd1d2.DeltaR(fJ2)) < 0.8 && fabs(Genpartd2.DeltaR(fJ2)) < 0.8 ){
                matching[1] = 5;
            }
            if( ( (fabs(Genpartd1d1.DeltaR(fJ2)) < 0.8)+(fabs(Genpartd1d2.DeltaR(fJ2)) < 0.8) )==1 && fabs(Genpartd2.DeltaR(fJ2)) < 0.8 ){
                matching[1] = 6;
            }
        }

    }
}

// 1 : matching with GenW
// 2 : matching with GenW, d1, d2
// 3 : matching with q,g
// 4 : matching with GenW d1 or d2
// 5 : matching with T, W, b, Wd1, Wd2
// 6 : matching with T, W, b, Wd1 or Wd2

std::array<int, 2> matching_array;
matching_array.at(0) = matching[0];
matching_array.at(1) = matching[1];
return matching_array;
'''

corr_dnn = '''
double dnn_wMjLL[19] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
double dnn_wMjLH[19] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
// LL
double dnn_corr_t2_wMjLL[18] = {1.00220259722,0.972082471336,0.956285915353,0.990345789614,0.978947300963,0.968049884931,0.995938115258,0.997333988826,1.03234025044,1.06018610614,1.06899448576,1.08514810247,1.23598621477,1.2502728724,1.30890928382,1.36304009108,1.38648970433,1.57352692135};
double dnn_corr_w_wMjLL[18] = {2.19914393291,2.03010441312,1.93473183745,1.81640348944,1.77027121281,1.66297099829,1.66342538332,1.36771111349,1.29514838164,1.16030381188,1.12637927884,1.02414709047,1.0046802008,0.955555341348,0.913972166178,0.816438221352,0.722388973257,0.583196347705};
double dnn_corr_g_wMjLL[18] = {0.573572134366,0.750227253152,0.874032113502,0.955071613117,1.02550960091,1.09549196544,1.15741760191,1.20831182514,1.25759903805,1.32870923916,1.35292783625,1.39982068505,1.38515812299,1.36379236233,1.29360026073,1.19993457302,1.09487314411,0.882450069856};

// LH
double dnn_corr_t2_wMjLH[18] = {0.956319424007,0.945284829286,0.931714248013,0.923526720708,0.994513733067,0.990917944939,0.989894236154,0.987249087179,1.11494668957,1.12852761348,0.99707201956,1.02157140796,1.01850430772,1.0225539138,1.0132037891,1.0779933365,1.08765838427,1.23099793976};
double dnn_corr_w_wMjLH[18] = {2.91915410287,2.71861787262,2.47199411528,2.32319881267,2.26957109897,2.2093226998,2.19217018322,1.53342581544,1.39335128866,1.022732901,0.995792746319,0.887309363694,0.880051657502,0.858984562377,0.853719906185,0.783933012948,0.758302898648,0.573342286031};
double dnn_corr_g_wMjLH[18] = {0.641635054671,0.788366997801,0.887441215083,0.957819854548,1.00939518637,1.06618708449,1.10690247535,1.15329413544,1.21277821534,1.28075774299,1.30476628105,1.33175221504,1.33875494495,1.34685089237,1.327900675,1.30379234562,1.26054335202,1.25505411938};

float Mj_1      = Lep1fatJet2_FatJet_msoftdrop;
float Mj_2      = Lep1fatJet2_FatJet_msoftdrop_2;
float PTj_1     = Lep1fatJet2_FatJet_pt;
float PTj_2     = Lep1fatJet2_FatJet_pt_2;
float dnnW_MD_1 = Lep1fatJet2_FatJet_deepTagMD_WvsQCD;
float dnnW_MD_2 = Lep1fatJet2_FatJet_deepTagMD_WvsQCD_2;

int matching_W = 2;
int matching_Wd = 4;
int matching_qg = 3;
int matching_T2 = 6;

float corr_deep[2] = {1,1};

// =================== correct genW, qg ================
// =================== correct genW, qg ================
// =================== correct genW, qg ================
// LL
if( Mj_1>60 && Mj_1<120 && PTj_1<400 && PTj_1>200 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_1>dnn_wMjLL[i] && dnnW_MD_1<dnn_wMjLL[i+1]){
            if( matching[0] == matching_T2 ){
                corr_deep[0] = corr_deep[0]*dnn_corr_t2_wMjLL[i];
            }
            if( matching[0] == matching_W ){
                corr_deep[0] = corr_deep[0]*dnn_corr_w_wMjLL[i];
            }
            if( matching[0] == matching_qg ){
                corr_deep[0] = corr_deep[0]*dnn_corr_g_wMjLL[i];
            }
        }
    }
}

if( Mj_2>60 && Mj_2<120 && PTj_2<400 && PTj_2>200 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_2>dnn_wMjLL[i] && dnnW_MD_2<dnn_wMjLL[i+1]){
            if( matching[1] == matching_T2 ){
                corr_deep[0] = corr_deep[0]*dnn_corr_t2_wMjLL[i];
            }
            if( matching[1] == matching_W ){
                corr_deep[0] = corr_deep[0]*dnn_corr_w_wMjLL[i];
            }
            if( matching[1] == matching_qg ){
                corr_deep[0] = corr_deep[0]*dnn_corr_g_wMjLL[i];
            }
        }
    }
}

// LH
if( Mj_1>60 && Mj_1<120 && PTj_1>400 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_1>dnn_wMjLH[i] && dnnW_MD_1<dnn_wMjLH[i+1]){
            if( matching[0] == matching_T2 ){
                corr_deep[0] = corr_deep[0]*dnn_corr_t2_wMjLH[i];
            }
            if( matching[0] == matching_W ){
                corr_deep[0] = corr_deep[0]*dnn_corr_w_wMjLH[i];
            }
            if( matching[0] == matching_qg ){
                corr_deep[0] = corr_deep[0]*dnn_corr_g_wMjLH[i];
            }
        }
    }
}

if( Mj_2>60 && Mj_2<120 && PTj_2>400 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_2>dnn_wMjLH[i] && dnnW_MD_2<dnn_wMjLH[i+1]){
            if( matching[1] == matching_T2 ){
                corr_deep[0] = corr_deep[0]*dnn_corr_t2_wMjLH[i];
            }
            if( matching[1] == matching_W ){
                corr_deep[0] = corr_deep[0]*dnn_corr_w_wMjLH[i];
            }
            if( matching[1] == matching_qg ){
                corr_deep[0] = corr_deep[0]*dnn_corr_g_wMjLH[i];
            }
        }
    }
}

// =================== correct genW, qg ================
// =================== correct genW, qg ================
// =================== correct genW, qg ================


// =================== correct genW, qg, genW daughter =================== 
// =================== correct genW, qg, genW daughter =================== 
// =================== correct genW, qg, genW daughter =================== 
// LL
if( Mj_1>60 && Mj_1<120 && PTj_1<400 && PTj_1>200 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_1>dnn_wMjLL[i] && dnnW_MD_1<dnn_wMjLL[i+1]){
            if( matching[0] == matching_T2 ){
                corr_deep[1] = corr_deep[1]*dnn_corr_t2_wMjLL[i];
            }
            if( matching[0] == matching_W ){
                corr_deep[1] = corr_deep[1]*dnn_corr_w_wMjLL[i];
            }
            if( matching[0] == matching_qg ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLL[i];
            }
            if( matching[0] == matching_Wd ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLL[i];
            }
        }
    }
}

if( Mj_2>60 && Mj_2<120 && PTj_2<400 && PTj_2>200 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_2>dnn_wMjLL[i] && dnnW_MD_2<dnn_wMjLL[i+1]){
            if( matching[1] == matching_T2 ){
                corr_deep[1] = corr_deep[1]*dnn_corr_t2_wMjLL[i];
            }
            if( matching[1] == matching_W ){
                corr_deep[1] = corr_deep[1]*dnn_corr_w_wMjLL[i];
            }
            if( matching[1] == matching_qg ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLL[i];
            }
            if( matching[1] == matching_Wd ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLL[i];
            }
        }
    }
}

// LH
if( Mj_1>60 && Mj_1<120 && PTj_1>400 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_1>dnn_wMjLH[i] && dnnW_MD_1<dnn_wMjLH[i+1]){
            if( matching[0] == matching_T2 ){
                corr_deep[1] = corr_deep[1]*dnn_corr_t2_wMjLH[i];
            }
            if( matching[0] == matching_W ){
                corr_deep[1] = corr_deep[1]*dnn_corr_w_wMjLH[i];
            }
            if( matching[0] == matching_qg ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLH[i];
            }
            if( matching[0] == matching_Wd ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLH[i];
            }
        }
    }
}

if( Mj_2>60 && Mj_2<120 && PTj_2>400 ){
    for(int i=0;i<18;i++ ){
        if( dnnW_MD_2>dnn_wMjLH[i] && dnnW_MD_2<dnn_wMjLH[i+1]){
            if( matching[1] == matching_T2 ){
                corr_deep[1] = corr_deep[1]*dnn_corr_t2_wMjLH[i];
            }
            if( matching[1] == matching_W ){
                corr_deep[1] = corr_deep[1]*dnn_corr_w_wMjLH[i];
            }
            if( matching[1] == matching_qg ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLH[i];
            }
            if( matching[1] == matching_Wd ){
                corr_deep[1] = corr_deep[1]*dnn_corr_g_wMjLH[i];
            }
        }
    }
}

// =================== correct genW, qg, genW daughter =================== 
// =================== correct genW, qg, genW daughter =================== 
// =================== correct genW, qg, genW daughter =================== 


std::array<float, 2> corr_dnn;
for(size_t ik=0; ik<corr_dnn.size();ik++){
    corr_dnn.at(ik) = corr_deep[ik];
}

return corr_dnn;

// corr_dnn[0] : correct genW, qg
// corr_dnn[1] : correct genW, qg, single genW daughter
'''

