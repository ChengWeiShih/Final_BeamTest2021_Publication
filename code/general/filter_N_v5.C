struct hit_str{
    int adc;
    int ampl;
    int chip_id;
    int fpga_id;
    int module;
    int chan_id;
    int fem_id;
    int bco;
    int bco_full;
    int event;
};

void filter_N_v5 (TString file_name) {

	// note : ===================================option===========================================
	// TString file_name = "BeamData_20211210-2043_0";
	// TString file_name = "run71_no_clone";
	// TString path_in = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/Run52/";
	
	TString path_in = "/data4/chengwei/Geant4/INTT_simulation/G4/for_CW/Final_BeamTest2021_Publication/data/DAC_Scan/";
	// note : ===================================option===========================================
	
	TFile * file_in = new TFile(path_in+file_name+".root", "read");
	TTree * tree_both_in = (TTree*)file_in -> Get("tree_both");

    vector<int> *vec_camacadc_bo = 0;
	vector<int> *vec_camactdc_bo = 0;
	bool INTT_event_bo;

	vector<int> * vec_adc=0;
	vector<int> * vec_ampl=0;
	vector<int> * vec_chip_id=0;
	vector<int> * vec_fpga_id=0;
	vector<int> * vec_module=0;
	vector<int> * vec_chan_id=0;
	vector<int> * vec_fem_id=0;
	vector<int> * vec_bco=0;
	vector<int> * vec_bco_full=0;
	vector<int> * vec_event=0;
	

	tree_both_in -> SetBranchAddress("camac_adc", &vec_camacadc_bo);
	tree_both_in -> SetBranchAddress("camac_tdc", &vec_camactdc_bo);
	tree_both_in -> SetBranchAddress("INTT_event", &INTT_event_bo);

	tree_both_in -> SetBranchAddress("adc", &vec_adc);
	tree_both_in -> SetBranchAddress("ampl", &vec_ampl);
	tree_both_in -> SetBranchAddress("chip_id", &vec_chip_id);
	tree_both_in -> SetBranchAddress("fpga_id", &vec_fpga_id);
	tree_both_in -> SetBranchAddress("module", &vec_module);
	tree_both_in -> SetBranchAddress("chan_id", &vec_chan_id);
	tree_both_in -> SetBranchAddress("fem_id", &vec_fem_id);
	tree_both_in -> SetBranchAddress("bco", &vec_bco);
	tree_both_in -> SetBranchAddress("bco_full", &vec_bco_full);
	tree_both_in -> SetBranchAddress("event", &vec_event);

	// note : # of the event in "tree_both"
	long nEvent_both = tree_both_in -> GetEntriesFast();

	Printf("tree_both,	# of event in %s.root : [%li]\n",file_name.Data(),nEvent_both);

	// note : ======================preparation of recreate================================
	// note : the module ID of each layer
	int module_l0 =1;
	int module_l1 =6;
	int module_l2 =5;
	int module_l3 =100; //note : don't need to see this 


	vector<int> * camac_adc_out;
	vector<int> * camac_tdc_out;
	bool INTT_event_out;

	vector<int> adc_array;      adc_array.clear();
	vector<int> ampl_array;     ampl_array.clear();
	vector<int> chip_id_array;  chip_id_array.clear();
	vector<int> fpga_id_array;  fpga_id_array.clear();
	vector<int> module_array;   module_array.clear();
	vector<int> chan_id_array;  chan_id_array.clear();
	vector<int> fem_id_array;   fem_id_array.clear();
	vector<int>	bco_array;      bco_array.clear();
	vector<int> bco_full_array; bco_full_array.clear();
	vector<int> event_array;    event_array.clear();
	int vec_element;
	int DSE_element = 0; // note : the default is 0 that means "not a double saving event"
	int eID_element = 0; // note : the true event ID, from 0 ~ ....

	int parent_eID = 0;
	int sub_eID    = 0;

	TFile *file_output = new TFile(path_in+file_name+"_filter_all_v5.root", "RECREATE");
	TTree *tree_output3 = new TTree("tree_both", "tree_both");
	
	tree_output3->Branch("camac_adc"        ,&camac_adc_out);
	tree_output3->Branch("camac_tdc"        ,&camac_tdc_out);
	tree_output3->Branch("INTT_event"       ,&INTT_event_out);
	tree_output3->Branch("adc"              ,&adc_array);
	tree_output3->Branch("ampl"             ,&ampl_array);
	tree_output3->Branch("chip_id"          ,&chip_id_array);
	tree_output3->Branch("fpga_id"          ,&fpga_id_array);
	tree_output3->Branch("module"           ,&module_array);
	tree_output3->Branch("chan_id"          ,&chan_id_array);
	tree_output3->Branch("fem_id"           ,&fem_id_array);
	tree_output3->Branch("bco"              ,&bco_array);
	tree_output3->Branch("bco_full"         ,&bco_full_array);
	tree_output3->Branch("event"            ,&event_array);
	tree_output3->Branch("nele"             ,&vec_element);
	tree_output3->Branch("DSE"              ,&DSE_element); // note : Double Saving Event check 
	tree_output3->Branch("eID"              ,&eID_element); // note : the line to save the real eID for quick check

    std::map<int, std::vector<hit_str>> hit_seg1_map; hit_seg1_map.clear();
    std::vector<hit_str> empty_hit_vec; empty_hit_vec.clear();

    for (int i=0; i<nEvent_both; i++) {
        tree_both_in -> GetEntry(i);

        parent_eID = i;

        if (i%10==0){cout<<"we are working on tree_both filter : "<<i<<endl;}

        unsigned int nhit_evt = vec_adc->size();
        if (nhit_evt == 0) {continue;} // note: skip empty events

        camac_adc_out = vec_camacadc_bo;
		camac_tdc_out = vec_camactdc_bo;
		INTT_event_out = INTT_event_bo;

        hit_seg1_map.clear();
        int BcoBcoFull_diff_count = 0;

        int first_bco = vec_bco->at(0);
        int first_bco_full = vec_bco_full->at(0);

        hit_seg1_map.insert(
            std::make_pair(
                BcoBcoFull_diff_count,
                // Form("%d_%d_%d", BcoBcoFull_diff_count, first_bco, first_bco_full),
                empty_hit_vec
            )
        );

        for (int hit_i = 0; hit_i < nhit_evt; hit_i++) {
            hit_str hit;
            hit.adc = vec_adc->at(hit_i);
            hit.ampl = vec_ampl->at(hit_i);
            hit.chip_id = vec_chip_id->at(hit_i);
            hit.fpga_id = vec_fpga_id->at(hit_i);
            hit.module = vec_module->at(hit_i);
            hit.chan_id = vec_chan_id->at(hit_i);
            hit.fem_id = vec_fem_id->at(hit_i);
            hit.bco = vec_bco->at(hit_i);
            hit.bco_full = vec_bco_full->at(hit_i);
            hit.event = vec_event->at(hit_i);

            if (hit.bco != first_bco || hit.bco_full != first_bco_full) {
                BcoBcoFull_diff_count += 1;
                first_bco = hit.bco;
                first_bco_full = hit.bco_full;

                hit_seg1_map.insert(
                    std::make_pair(
                        BcoBcoFull_diff_count,
                        // Form("%d_%d_%d", BcoBcoFull_diff_count, hit.bco, hit.bco_full),
                        empty_hit_vec
                    )
                );
            }

            // if (hit_seg1_map.find(Form("%d_%d_%d", BcoBcoFull_diff_count, hit.bco, hit.bco_full)) != hit_seg1_map.end()) {
            //     hit_seg1_map[Form("%d_%d_%d", BcoBcoFull_diff_count, hit.bco, hit.bco_full)].push_back(hit);
            // }

            if (hit_seg1_map.find(BcoBcoFull_diff_count) != hit_seg1_map.end()) {
                hit_seg1_map[BcoBcoFull_diff_count].push_back(hit);
            }

            else {
                std::cout << "Error: Key not found in hit_seg1_map: " << BcoBcoFull_diff_count << std::endl;
                exit(1);
            }
        } // note : end of nhit

        sub_eID = 0;
        for (const auto& pair : hit_seg1_map) { // note : looping each sub event
            adc_array.clear();
            ampl_array.clear();
            chip_id_array.clear();
            fpga_id_array.clear();
            module_array.clear();
            chan_id_array.clear();
            fem_id_array.clear();
            bco_array.clear();
            bco_full_array.clear();
            event_array.clear();

            // const std::string& key = pair.first;
            const int& key = pair.first;
            const std::vector<hit_str>& hits = pair.second;
            std::map<std::string, hit_str> SubEvt_hit_map; SubEvt_hit_map.clear();

            if (hits.size() == 0) {
                std::cout<< "Warning: No hits found for key: " << key << std::endl;
                exit(1);
            }

            for (int hit_i = 0; hit_i < hits.size(); hit_i++) { // note : hits in each sub event

                if (hits[hit_i].module != module_l0 && hits[hit_i].module != module_l1 && hits[hit_i].module != module_l2) {continue;}
                if (hits[hit_i].chip_id < 1 || hits[hit_i].chip_id > 26) {continue;} // note : range : 1 ~ 26
                if (hits[hit_i].chan_id < 0 || hits[hit_i].chan_id > 127) {continue;} // note : range : 0 ~ 127
                if (hits[hit_i].adc < 0 || hits[hit_i].adc > 7) {continue;} // note : range : 0 ~ 7
                if (hits[hit_i].bco < 0 || hits[hit_i].bco > 127) {continue;} // note : range : 0 ~ 127
                if (hits[hit_i].bco_full < 0 || hits[hit_i].bco_full > 65535) {continue;} // note : range : 0 ~ 65535

                // note : the latter will overwrite the former if the key is the same
                SubEvt_hit_map[Form("%d_%d_%d_%d_%d_%d", 
                    hits[hit_i].module, 
                    hits[hit_i].chip_id, 
                    hits[hit_i].chan_id, 
                    hits[hit_i].adc, 
                    hits[hit_i].bco, 
                    hits[hit_i].bco_full
                )] = hits[hit_i];                
            } // note : end of hits in each sub event

            for (auto &pair : SubEvt_hit_map){
                hit_str& hit = pair.second;

                adc_array.push_back(hit.adc);
                ampl_array.push_back(hit.ampl);
                chip_id_array.push_back(hit.chip_id);
                fpga_id_array.push_back(hit.fpga_id);
                module_array.push_back(hit.module);
                chan_id_array.push_back(hit.chan_id);
                fem_id_array.push_back(hit.fem_id);
                bco_array.push_back(hit.bco);
                bco_full_array.push_back(hit.bco_full);
                event_array.push_back(hit.event);
            }

            vec_element = adc_array.size();
            DSE_element = 0; // note : the clone hits are removed by the map already
            eID_element = 1000*parent_eID+sub_eID;

            tree_output3->Fill();
            // if (vec_element != 0) {tree_output3->Fill();}

            sub_eID += 1;

            if (sub_eID > 999){
                std::cout<<"!!!!!!! eID re-formating error !!!!!!!!! "<<", sub_eID: "<<sub_eID<<", eID re-formating: "<<eID_element<<std::endl;
            }

        } // note : end of hit_seg1_map

    } // note : end of Event


    cout<<"test test "<<endl;
	tree_output3->Write("", TObject::kOverwrite);
	cout<<"test test "<<endl;
	file_output->Close();
	printf("filter done\n");
	file_in->Close();
}