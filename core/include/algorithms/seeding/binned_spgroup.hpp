/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

class binned_spgroup_iterator{
private:
    
}

class binned_spgroup{
    binned_spgroup(const spacepoint_container& sp_container,
		   const spacepoint_grid& sp_grid){
	m_sp_container = sp_container;
	m_sp_grid = sp_grid;
    }    
    
private:
    spacepoint_container m_sp_container;
    spacepoint_grid m_sp_grid;

};
