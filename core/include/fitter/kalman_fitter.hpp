/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

namespace traccc {

struct kalman_fitter_option{


};

template <typename propagator_t >    
class kalman_fitter {
public:
    using navigatition_surface = typename propagator_t::navigation_surface;

    /// Default constructor is deleted
    kalman_fitter() = delete;
    
    /// Constructor from arguments
    kalman_fitter(propagator_t propagator)
      : m_propagator(std::move(propagator)) {}


private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

    
};
    

} // namespace traccc    
