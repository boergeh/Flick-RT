#ifndef flick_single_layer_slab
#define flick_single_layer_slab

#include "../numeric/named_bounded_types.hpp"
#include "../material/henyey_greenstein.hpp"
#include "../transporter/ordinary_mc.hpp"

namespace flick {
namespace model {
  class single_layer_slab
  {
    thickness h_;
    zenith_angle theta_0_{0};
    albedo albedo_{0};
    semi_infinite_box geometry_;
    number n_packages_{1};
    std::shared_ptr<material::base> material_;
    receiver* transmitted_;
    receiver* reflected_;
    double relative_depth_{0};
    stokes stokes_{stokes::unpolarized()};
    std::shared_ptr<transporter::ordinary_mc> omc_;
  public:
    single_layer_slab(const thickness& h) : h_{h} {
    }
    void brighten_bottom(const albedo& a) {
      albedo_ = a;
    }
    void load_packages(const number& n) {
      n_packages_ = n;
    }
    void orient_source(const zenith_angle& za) {
      theta_0_ = za;
    }
    void initiate_source(const stokes& s) {
      stokes_ = s;
    }
    template <class Material, class... Args>
    void fill(Args... a) {
      material_ = std::make_shared<Material>(a...);
    }
    double hemispherical_reflectance(const unit_interval& relative_depth
				     = unit_interval{0})
    // See Wikipedia reflectance for definition
    {
      relative_depth_ = relative_depth();
      run();
      find_receivers();
      return reflected_->radiant_flux()/n_packages_();
    }
    double relative_radiance(const polar_angle& pa,
			     const azimuth_angle& aa,
			     const vertex_angle& acceptance_angle
			     = vertex_angle{1},
			     const unit_interval& relative_depth
			     = unit_interval{0})
    // Radiance in a given propagation direction relative to incoming
    // surface irradiance.
    {
      relative_depth_ = relative_depth();
      run();
      find_receivers();
      unit_vector direction{pa(),aa()};
      double L_r = reflected_->radiance(direction, acceptance_angle());
      double L_t = transmitted_->radiance(direction, acceptance_angle());
      return (L_r + L_t) / n_packages_();
    }
    double hemispherical_transmittance(const unit_interval& relative_depth
				       = unit_interval{1})
    // See Wikipedia transmittance for definition
    {
      relative_depth_ = relative_depth();
      run();
      find_receivers();
      return transmitted_->radiant_flux()/n_packages_();
    }
  private:
    double relative_skin_depth() {
      return geometry_.small_step()/h_()*2;
    }
    double sheet_relative_depth() {
      double epsilon = relative_skin_depth();
      return std::clamp<double>(relative_depth_,epsilon,1-epsilon);
    }
    void find_receivers() {
      transmitted_ = &omc_->inward_receiver("sheet");
      reflected_ = &omc_->outward_receiver("sheet");

      double epsilon = relative_skin_depth();
      if (relative_depth_ <= epsilon) {
	reflected_ = &omc_->outward_receiver("surface");
      }
      else if (relative_depth_ >= 1-epsilon) {
	transmitted_ = &omc_->inward_receiver("bottom");
      }
    }
    void build_geometry() {
      semi_infinite_box surface;
      semi_infinite_box sheet;
      semi_infinite_box bottom;
      geometry_.name("geometry");
      surface.name("surface");
      sheet.name("sheet");
      bottom.name("bottom");
      surface().inward_receiver().activate();
      surface().outward_receiver().activate();
      surface().fill(material_);
      sheet().inward_receiver().activate();
      sheet().outward_receiver().activate();
      sheet().fill(material_);
      bottom().inward_receiver().activate();
      bottom().coat<coating::grey_lambert>(albedo_(),1-albedo_());
      geometry_.move_by({0,0,h_()+1});
      surface.move_by({0,0,h_()});
      sheet.move_by({0,0,h_()*(1-sheet_relative_depth())});
      sheet.insert(bottom);
      surface.insert(sheet);
      geometry_.clear();
      geometry_.insert(surface);
    }
    void run() {
      emitter emitter{{0,0,h_()+0.5},stokes_,n_packages_()};
      unit_vector direction{constants::pi-theta_0_(),0}; 
      emitter.set_direction<unidirectional>(direction);
      build_geometry();
      double g = material_->asymmetry_factor();
      omc_ = std::make_shared<transporter::ordinary_mc>(geometry_);
      omc_->transport_radiation(emitter,"geometry",g);
    }
  };
}
}

#endif
