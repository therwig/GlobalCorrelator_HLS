#include "tmux_create_test.h"

int mp7DataLength = 2*(NTRACK+NCALO+NEMCALO+NMU);
int objDataLength[4] = {2*NTRACK, 2*(NEMCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK), 2*(NEMCALO+NCALO+NTRACK+NMU)};
int link_max[4] = {NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_CALO+NLINKS_PER_EMCALO+NLINKS_PER_MU};
int link_min[4] = {0, NLINKS_PER_TRACK, NLINKS_PER_TRACK+NLINKS_PER_EMCALO, NLINKS_PER_TRACK+NLINKS_PER_EMCALO+NLINKS_PER_CALO};
unsigned int theEtaRegion = 0;
unsigned int thePhiRegion = 0;

int main() {

    if (theEtaRegion>=NETA_TMUX) theEtaRegion = NETA_TMUX-1;
    if (thePhiRegion>=NPHI_TMUX) thePhiRegion = NPHI_TMUX-1;

    // input format: could be random or coming from simulation
    //RandomPFInputs inputs(37); // 37 is a good random number
    //DiscretePFInputs inputs("/home/therwig/sandbox/GlobalCorrelator_HLS/data/regions_TTbar_PU140.dump");
    //DiscretePFInputs inputs("barrel_sectors_1x1_TTbar_PU140.dump");
    DiscretePFInputs inputs("/home/therwig/sandbox/GlobalCorrelator_HLS/data/barrel_sectors_1x1_TTbar_PU200.dump");
    //DiscretePFInputs inputs("dummy.dump");
    
    // input TP objects
    HadCaloObj calo[NCALO_TMUX]; EmCaloObj emcalo[NEMCALO_TMUX]; TkObj track[NTRACK_TMUX]; z0_t hwZPV;
    HadCaloObj calo_subem[NCALO_TMUX], calo_subem_ref[NCALO_TMUX]; 
    MuObj mu[NMU_TMUX];

    //printf("NTRACK = %i, NEMCALO = %i, NCALO = %i, NMU = %i, MP7_NCHANN = %i \n", NTRACK, NEMCALO, NCALO, NMU, MP7_NCHANN);

    //std::cout<<mp7DataLength<<std::endl;
    const int listLength = NFRAMES_APX_GEN0*((NTEST*TMUX_OUT)+(TMUX_IN-TMUX_OUT));
    // 3 * (6*18 + 36-18) = 3*126 = 378
    //std::cout<<listLength<<std::endl;
    std::string datawords [NLINKS_APX_GEN0][listLength];
    for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
        for (int ib = 0; ib < listLength; ib++){
            datawords[ia][ib] = "0x0000000000000000";
        }
    }

    // -----------------------------------------
    // run multiple tests

    int link_off = 0;
    int offset = 0;
    for (int test = 0; test < NTEST; ++test) {
        

        // initialize TP objects
        for (int i = 0; i < NTRACK_TMUX; ++i) { // 15 * 36 * 1 * 2 = 1080
            track[i].hwPt = 0; track[i].hwPtErr = 0; track[i].hwEta = 0; track[i].hwPhi = 0; track[i].hwZ0 = 0; 
        }
        for (int i = 0; i < NCALO_TMUX; ++i) {
            calo[i].hwPt = 0; calo[i].hwEmPt = 0; calo[i].hwEta = 0; calo[i].hwPhi = 0; calo[i].hwIsEM = 0; 
        }
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            emcalo[i].hwPt = 0; emcalo[i].hwPtErr = 0;  emcalo[i].hwEta = 0; emcalo[i].hwPhi = 0;
        }
        for (int i = 0; i < NMU_TMUX; ++i) {
            mu[i].hwPt = 0; mu[i].hwPtErr = 0; mu[i].hwEta = 0; mu[i].hwPhi = 0;
        }

        // get the inputs from the input object
        if (!inputs.nextRegion_tmux(calo, emcalo, track, mu, hwZPV)) break;

        /*for (int i = 0; i < NTRACK_TMUX; ++i) {
            std::cout<<track[i].hwPt<<"\t "<<track[i].hwEta<<"\t "<<track[i].hwPhi<<std::endl;
        }
        for (int i = 0; i < NCALO_TMUX; ++i) {
            std::cout<<calo[i].hwPt<<"\t "<<calo[i].hwEta<<"\t "<<calo[i].hwPhi<<std::endl;
        }
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            std::cout<<emcalo[i].hwPt<<"\t "<<emcalo[i].hwEta<<"\t "<<emcalo[i].hwPhi<<std::endl;
        }
        for (int i = 0; i < NMU_TMUX; ++i) {
            std::cout<<mu[i].hwPt<<"\t "<<mu[i].hwEta<<"\t "<<mu[i].hwPhi<<std::endl;
        }*/

        //VtxObj curvtx;    
        //simple_vtx_ref(track,&curvtx);
        //printf("Vertex Z   %i\n",(int)(curvtx.hwZ0));

        unsigned int ie = theEtaRegion;
        unsigned int ip = thePhiRegion;
        //std::cout<<"ie"<<ie<<" ip"<<ip<<std::endl;
        HadCaloObj calo_temp[TMUX_OUT][NCALO]; EmCaloObj emcalo_temp[TMUX_OUT][NEMCALO]; TkObj track_temp[TMUX_OUT][NTRACK]; MuObj mu_temp[TMUX_OUT][NMU];
        for (int ir = 0; ir < TMUX_OUT; ir++) {
            // initialize temp objects
            for (int i = 0; i < NTRACK; ++i) {
                track_temp[ir][i].hwPt = 0; track_temp[ir][i].hwPtErr = 0; track_temp[ir][i].hwEta = 0; track_temp[ir][i].hwPhi = 0; track_temp[ir][i].hwZ0 = 0; 
            }
            for (int i = 0; i < NCALO; ++i) {
                calo_temp[ir][i].hwPt = 0; calo_temp[ir][i].hwEmPt = 0; calo_temp[ir][i].hwEta = 0; calo_temp[ir][i].hwPhi = 0; calo_temp[ir][i].hwIsEM = 0; 
            }
            for (int i = 0; i < NEMCALO; ++i) {
                emcalo_temp[ir][i].hwPt = 0; emcalo_temp[ir][i].hwPtErr = 0;  emcalo_temp[ir][i].hwEta = 0; emcalo_temp[ir][i].hwPhi = 0;
            }
            for (int i = 0; i < NMU; ++i) {
                mu_temp[ir][i].hwPt = 0; mu_temp[ir][i].hwPtErr = 0; mu_temp[ir][i].hwEta = 0; mu_temp[ir][i].hwPhi = 0;
            }
        }
        // fill temp containers
        int etalo = -MAXETA_INT+int(float(2*MAXETA_INT*ie)/float(NETA_TMUX))-ETA_BUFFER;
        int etahi = -MAXETA_INT+int(float(2*MAXETA_INT*(ie+1))/float(NETA_TMUX))+ETA_BUFFER;
        int philo = -MAXPHI_INT+int(float(2*MAXPHI_INT*ip)/float(NPHI_TMUX))-PHI_BUFFER;
        int phihi = -MAXPHI_INT+int(float(2*MAXPHI_INT*(ip+1))/float(NPHI_TMUX))+PHI_BUFFER;
        //std::cout<<etalo<<" "<<etahi<<" "<<philo<<" "<<phihi<<" "<<std::endl;

        int i_temp = 0;
        int ireg = 0;
        int ntracks = 0;
        int ncalos = 0;
        int nemcalos = 0;
        int nmus = 0;
        for (int i = 0; i < NTRACK_TMUX; ++i) {
            if (int(track[i].hwEta) < etalo or int(track[i].hwEta) > etahi) continue;
            if (int(track[i].hwPhi) < philo or int(track[i].hwPhi) > phihi) continue;
            if (int(track[i].hwPt) == 0) continue;
            track_temp[ireg][i_temp] = track[i];
            ntracks++;
            ireg++;
            if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
            if (i_temp == NTRACK) {break;}
        }
        i_temp = 0;
        ireg = 0;
        for (int i = 0; i < NCALO_TMUX; ++i) {
            if (int(calo[i].hwEta) < etalo or int(calo[i].hwEta) > etahi) continue;
            if (int(calo[i].hwPhi) < philo or int(calo[i].hwPhi) > phihi) continue;
            if (int(calo[i].hwPt) == 0) continue;
            calo_temp[ireg][i_temp] = calo[i];
            ncalos++;
            ireg++;
            if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
            if (i_temp == NCALO) {break;}
        }
        i_temp = 0;
        ireg = 0;
        for (int i = 0; i < NEMCALO_TMUX; ++i) {
            if (int(emcalo[i].hwEta) < etalo or int(emcalo[i].hwEta) > etahi) continue;
            if (int(emcalo[i].hwPhi) < philo or int(emcalo[i].hwPhi) > phihi) continue;
            if (int(emcalo[i].hwPt) == 0) continue;
            emcalo_temp[ireg][i_temp] = emcalo[i];
            nemcalos++;
            ireg++;
            if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
            if (i_temp == NEMCALO) {break;}
        }
        i_temp = 0;
        ireg = 0;
        /*for (int i = 0; i < NMU_TMUX; ++i) {
            if (int(mu[i].hwEta) < etalo or int(mu[i].hwEta) > etahi) continue;
            if (int(mu[i].hwPhi) < philo or int(mu[i].hwPhi) > phihi) continue;
            if (int(mu[i].hwPt) == 0) continue;
            mu_temp[ireg][i_temp] = mu[i];
            nmus++;
            ireg++;
            if (ireg == TMUX_OUT) {i_temp++; ireg=0;}
            if (i_temp == NMU) {break;}
        }*/
        //FIXME skipping muons for now, makes link alignment with firmware work out

        std::cout<<"\n\ntrack  = "<<ntracks<<std::endl;

        // std::cout<<"Totals for event "<<ir<<":"<<std::endl;
        // std::cout<<"\ttrack  = "<<ntracks<<std::endl;
        // std::cout<<"\tcalo   = "<<ncalos<<std::endl;
        // std::cout<<"\temcalo = "<<nemcalos<<std::endl;
        // std::cout<<"\tmu     = "<<nmus<<std::endl;

        for (int ir = 0; ir < TMUX_OUT; ir++) {

            MP7DataWord data_in[MP7_NCHANN];
            mp7wrapped_pack_in_reorder(emcalo_temp[ir], calo_temp[ir], track_temp[ir], mu_temp[ir], data_in);
            //for (unsigned int in = 0; in < MP7_NCHANN; in++){
            //    printf("data_in[%i] = %i \n", in, (int) data_in[in]);
            //}

	    // per-test offset
            offset = NFRAMES_APX_GEN0*TMUX_OUT*test;

            float tot_perc[4];
            int link_start[4];
            int add_off[4];

	    // which of the 10 links for tracks are to be used for this time-slice
	    // expressed as a 'percentage' and then an integer 
            tot_perc[0] = float(ir*NLINKS_PER_TRACK)/float(TMUX_OUT); // (1..18)/18 * 10
            tot_perc[1] = float(ir*NLINKS_PER_EMCALO)/float(TMUX_OUT);
            tot_perc[2] = float(ir*NLINKS_PER_CALO)/float(TMUX_OUT);
            tot_perc[3] = float(ir*NLINKS_PER_MU)/float(TMUX_OUT);
            link_start[0] = int(tot_perc[0]);
            link_start[1] = int(tot_perc[1]);
            link_start[2] = int(tot_perc[2]);
            link_start[3] = int(tot_perc[3]);
	    //offset (in links*frames) as a function of the time slice considered
	    // 3*36*(link pct) mod 3*36
            add_off[0] = int(float(NFRAMES_APX_GEN0*TMUX_IN)*tot_perc[0])%(NFRAMES_APX_GEN0*TMUX_IN);
            add_off[1] = int(float(NFRAMES_APX_GEN0*TMUX_IN)*tot_perc[1])%(NFRAMES_APX_GEN0*TMUX_IN);
            add_off[2] = int(float(NFRAMES_APX_GEN0*TMUX_IN)*tot_perc[2])%(NFRAMES_APX_GEN0*TMUX_IN);
            add_off[3] = int(float(NFRAMES_APX_GEN0*TMUX_IN)*tot_perc[3])%(NFRAMES_APX_GEN0*TMUX_IN);
	    /* 
	       First few entries in each for track case (nlink=10):
	       link fractions are (0,0.55,1.11,1.66,2.22,...) 
	       rounded to integers (which ACTUAL linkes) (0,0,1,1,2,...) 
	       and leftover bits are offset by frames, indexed by (0,60,12,72,24,...)
	    */

            int id = 0;
            unsigned int link_type = 0; //0=track, 1=emcalo, 2=calo, 3=mu

	    //loop over the 10+10+10+2 = 32 links per region
            for (int link_ctr = 0; link_ctr < NLINKS_PER_REG; link_ctr++) {

                if      (link_ctr < link_max[0]) link_type = 0;
                else if (link_ctr < link_max[1]) link_type = 1;
                else if (link_ctr < link_max[2]) link_type = 2;
                else if (link_ctr < link_max[3]) link_type = 3;

                //std::cout<<"Region "<<ir<<" Link "<<link_ctr<<" : "<<tot_perc[link_type]<<" "<<link_start[link_type]<<" "<<add_off[link_type]<<std::endl;

		// (**) advance to the particular link to be considered for this time slice: link start (0,10,20,...) + tmux offset (1..18)/18*10
                if (link_ctr < link_min[link_type]+link_start[link_type]) {continue;}

                std::stringstream stream1;
                int index = add_off[link_type]; // the "overlap mod(frame*link)" offset
                if (index==0) { // only need the vertex for the first entry on this link
                    stream1 << "0x";
                    stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]); // TODO what is this data?
                    stream1 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00"; // add z0 to the second string
		    // full 64b word written (2 * 8 * 3=log_2(8))
                    datawords[link_off+link_ctr][offset] = stream1.str();
		    // link off = (0,..)
		    // link_ctr see (**)
		    // offset is per-TEST (frames*TMUX_OUT)
                    id++;
                    index++;
                }
                else {
		    // TODO: apparently remove the 6 least significant bits and replace with vertex ??
                    stream1 << datawords[link_off+link_ctr][offset].substr(0,10);
                    stream1 << std::setfill('0') << std::setw(6) << std::hex << (((unsigned int)(hwZPV.range(9,0))) << 14) << "00";
                    datawords[link_off+link_ctr][offset] = stream1.str();
                }

                //std::cout<<"index="<<index<<" id"<<id<<std::endl;
        
                // put the data on the link number = link_ctr;
		// consider the rest of the frames from the offset 'add_off' up to the limit of 3*36=108
                while (index < NFRAMES_APX_GEN0*TMUX_IN) {
                    stream1.str("");
                    stream1 << "0x";
                    if (index%NFRAMES_APX_GEN0==0) { // pad every third frame with 32b data + 32b zeros?
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]) << "00000000";
                        id+=1;
                    }
                    else if (id == objDataLength[link_type]-1) { // if is already at the max number of objects that may be stored
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]) << "00000000";
                        id+=1;
                    }
                    else {
		        //write two 32b words
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id+1]);
                        stream1 << std::setfill('0') << std::setw(8) << std::hex << (unsigned int) (data_in[id]);
                        id+=2;
                    }
                    datawords[link_off+link_ctr][offset+index] = stream1.str();
                    index++;
		    // write to link "link_off+link_ctr" and "offset+index" out of "listLength"=378
		    //   here, sequentially deeper into the second index

                    // std::cout << stream1.str() << std::endl;
                    // std::cout << datawords[link_ctr][offset+(id/2)] << std::endl;
                    //printf("test = %i, link ctr = %i, clock = %i, offset = %i \n", test, link_ctr, offset+((id+1)/2), offset);
                    //std::cout<<"link_ctr = "<<link_ctr<<" link_off = "<<link_off<<" offset = "<<offset<<" index = "<<index;
                    //std::cout<<"  "<<datawords[link_off+link_ctr][offset+index-1]<<std::endl;
                    if (id >= objDataLength[link_type]) {
                        break;
                    }
                }
                if (id >= objDataLength[link_type]) link_ctr = link_max[link_type]-1;
		//advance to the next link (corresponding to a different object) if the MAX is reached
		//  (for loop will increment link immediately following this and will advance to the next object type)
            }
            
            //std::cout<<"-----------"<<std::endl;
        }
        link_off += NLINKS_PER_REG;
        if (link_off>= NLINKS_PER_REG*TMUX_IN/TMUX_OUT) link_off = 0;
        //std::cout<<"\t"<<test<<std::endl;
    
    }

    for (int ib = 0; ib < listLength; ib++){
        // std::cout << ib << " ";
        std::cout << "0x" << std::setfill('0') << std::setw(4) << std::hex << ib << "   " <<std::dec;
        for (int ia = 0; ia < NLINKS_APX_GEN0; ia++){
            //datawords[ia][ib] = "0x0000000000000000";
            std::cout << datawords[ia][ib] << "    ";
        }
        std::cout << std::endl;
    }

    return 0;
}
