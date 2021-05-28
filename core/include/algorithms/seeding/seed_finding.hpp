/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */


#pragma once

#include <edm/internal_spacepoint.hpp>
#include <edm/seed.hpp>
#include <algorithms/seeding/detail/seeding_config.hpp>
#include <algorithms/seeding/doublet_finding.hpp>
#include <algorithms/seeding/triplet_finding.hpp>
#include <algorithms/seeding/seed_filtering.hpp>
#include <iostream>

namespace traccc{
    
struct seed_finding{

seed_finding(vecmem::memory_resource& resource,
	     seedfinder_config& config,
	     const host_internal_spacepoint_container& isp_container,
	     experiment_cuts* exp_cuts = nullptr):
    m_doublet_finding(resource, config, isp_container),
    m_triplet_finding(config, isp_container),
    m_seed_filtering(isp_container, exp_cuts),
    m_isp_container(isp_container)
    {}
    
host_seed_collection operator()(){
    host_seed_collection seed_collection;
    this->operator()(seed_collection);
    
    return seed_collection;
}
        
void operator()(host_seed_collection& seeds){
    
    // iterate over grid bins
    for (size_t i=0; i<m_isp_container.headers.size(); ++i){
	
	auto& bin_information = m_isp_container.headers[i];
	auto& spM_collection = m_isp_container.items[i];
	
	/// iterate over middle spacepoints
	for (size_t j=0; j<spM_collection.size(); ++j){
	    
	    sp_location spM_location({i,j});
	   
	    // doublet search	    
	    auto doublets_mid_bot = m_doublet_finding(bin_information, spM_location, true);
	    if (doublets_mid_bot.empty()) continue;
	    
	    auto doublets_mid_top = m_doublet_finding(bin_information, spM_location, false);
	    if (doublets_mid_top.empty()) continue;

	    //std::cout << doublets_mid_bot.size() << "  " << doublets_mid_top.size() << std::endl;
	    
	    //host_triplet_collection triplets_per_spM;
	    std::vector<triplet> triplets_per_spM;
	    
	    for (auto mid_bot: doublets_mid_bot){
		std::vector<triplet> triplets = m_triplet_finding(mid_bot, doublets_mid_top);//, triplets_per_spM);
		triplets_per_spM.insert(std::end(triplets_per_spM), triplets.begin(), triplets.end());
	    }

	    //std::cout << triplets_per_spM.size() << std::endl;
	    
	    m_seed_filtering(triplets_per_spM, seeds);	    	    	    
	}	
    }
}
    
private:
    const host_internal_spacepoint_container& m_isp_container;
    doublet_finding m_doublet_finding;
    triplet_finding m_triplet_finding;
    seed_filtering m_seed_filtering;
};
} // namespace traccc
