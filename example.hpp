
// Code generated by stanc f556d0df
#include <stan/model/model_header.hpp>
namespace example_model_namespace {

template <typename T, typename S>
std::vector<T> resize_to_match__(std::vector<T>& dst, const std::vector<S>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, -1>
resize_to_match__(Eigen::Matrix<T, -1, -1>& dst, const Eigen::Matrix<T, -1, -1>& src) {
  dst.resize(src.rows(), src.cols());
  return dst;
}

template <typename T>
Eigen::Matrix<T, 1, -1>
resize_to_match__(Eigen::Matrix<T, 1, -1>& dst, const Eigen::Matrix<T, 1, -1>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, 1>
resize_to_match__(Eigen::Matrix<T, -1, 1>& dst, const Eigen::Matrix<T, -1, 1>& src) {
  dst.resize(src.size());
  return dst;
}
std::vector<double> to_doubles__(std::initializer_list<double> x) {
  return x;
}

std::vector<stan::math::var> to_vars__(std::initializer_list<stan::math::var> x) {
  return x;
}

inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}

inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}


using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math; 

static int current_statement__ = 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 8, column 2 to column 37)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 9, column 2 to column 36)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 10, column 2 to column 32)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 11, column 2 to column 34)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 12, column 2 to column 33)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 16, column 2 to column 21)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 19, column 2 to line 20, column 31)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 22, column 2 to column 38)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 2, column 2 to column 18)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 3, column 2 to column 17)",
                                                      " (in '/content/drive/MyDrive/Github_rep/colab_python/example.stan', line 4, column 2 to column 17)"};



class example_model : public model_base_crtp<example_model> {

 private:
  int pos__;
  int N;
  Eigen::Matrix<double, -1, 1> wave;
  Eigen::Matrix<double, -1, 1> flux;
 
 public:
  ~example_model() { }
  
  std::string model_name() const { return "example_model"; }
  
  example_model(stan::io::var_context& context__,
                unsigned int random_seed__ = 0,
                std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    typedef double local_scalar_t__;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "example_model_namespace::example_model";
    (void) function__;  // suppress unused var warning
    
    try {
      
      pos__ = 1;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      
      current_statement__ = 9;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 10;
      validate_non_negative_index("wave", "N", N);
      context__.validate_dims("data initialization","wave","double",
          context__.to_vec(N));
      wave = Eigen::Matrix<double, -1, 1>(N);
      
      {
        std::vector<local_scalar_t__> wave_flat__;
        current_statement__ = 10;
        assign(wave_flat__, nil_index_list(), context__.vals_r("wave"),
          "assigning variable wave_flat__");
        current_statement__ = 10;
        pos__ = 1;
        current_statement__ = 10;
        for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 10;
          assign(wave, cons_list(index_uni(sym1__), nil_index_list()),
            wave_flat__[(pos__ - 1)], "assigning variable wave");
          current_statement__ = 10;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 11;
      validate_non_negative_index("flux", "N", N);
      context__.validate_dims("data initialization","flux","double",
          context__.to_vec(N));
      flux = Eigen::Matrix<double, -1, 1>(N);
      
      {
        std::vector<local_scalar_t__> flux_flat__;
        current_statement__ = 11;
        assign(flux_flat__, nil_index_list(), context__.vals_r("flux"),
          "assigning variable flux_flat__");
        current_statement__ = 11;
        pos__ = 1;
        current_statement__ = 11;
        for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 11;
          assign(flux, cons_list(index_uni(sym1__), nil_index_list()),
            flux_flat__[(pos__ - 1)], "assigning variable flux");
          current_statement__ = 11;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 9;
      current_statement__ = 9;
      check_greater_or_equal(function__, "N", N, 1);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename T__>
  T__ log_prob(std::vector<T__>& params_r__, std::vector<int>& params_i__,
               std::ostream* pstream__ = 0) const {
    typedef T__ local_scalar_t__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "example_model_namespace::log_prob";
(void) function__;  // suppress unused var warning

    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    
    try {
      local_scalar_t__ cont;
      
      current_statement__ = 1;
      cont = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        cont = stan::math::lub_constrain(cont, -1000, 1000, lp__);
      } else {
        current_statement__ = 1;
        cont = stan::math::lub_constrain(cont, -1000, 1000);
      }
      local_scalar_t__ slope;
      
      current_statement__ = 2;
      slope = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        slope = stan::math::lub_constrain(slope, -100, 100, lp__);
      } else {
        current_statement__ = 2;
        slope = stan::math::lub_constrain(slope, -100, 100);
      }
      local_scalar_t__ amp;
      
      current_statement__ = 3;
      amp = in__.scalar();
      current_statement__ = 3;
      if (jacobian__) {
        current_statement__ = 3;
        amp = stan::math::lub_constrain(amp, 0, 1000, lp__);
      } else {
        current_statement__ = 3;
        amp = stan::math::lub_constrain(amp, 0, 1000);
      }
      local_scalar_t__ center;
      
      current_statement__ = 4;
      center = in__.scalar();
      current_statement__ = 4;
      if (jacobian__) {
        current_statement__ = 4;
        center = stan::math::lub_constrain(center, 0, 10, lp__);
      } else {
        current_statement__ = 4;
        center = stan::math::lub_constrain(center, 0, 10);
      }
      local_scalar_t__ width;
      
      current_statement__ = 5;
      width = in__.scalar();
      current_statement__ = 5;
      if (jacobian__) {
        current_statement__ = 5;
        width = stan::math::lub_constrain(width, 0, 10, lp__);
      } else {
        current_statement__ = 5;
        width = stan::math::lub_constrain(width, 0, 10);
      }
      {
        current_statement__ = 6;
        validate_non_negative_index("mod_flux", "N", N);
        Eigen::Matrix<local_scalar_t__, -1, 1> mod_flux;
        mod_flux = Eigen::Matrix<local_scalar_t__, -1, 1>(N);
        
        current_statement__ = 6;
        for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 6;
          assign(mod_flux, cons_list(index_uni(sym1__), nil_index_list()),
            std::numeric_limits<double>::quiet_NaN(),
            "assigning variable mod_flux");}
        current_statement__ = 7;
        assign(mod_flux, nil_index_list(),
          add(
            add(
              multiply(amp,
                stan::math::exp(
                  divide(multiply(-0.5, square(subtract(center, wave))),
                    square(width)))), cont), multiply(slope, wave)),
          "assigning variable mod_flux");
        current_statement__ = 8;
        lp_accum__.add(
          normal_log<propto__>(flux, mod_flux, stan::math::sqrt(flux)));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob() 
    
  template <typename RNG>
  void write_array(RNG& base_rng__, std::vector<double>& params_r__,
                   std::vector<int>& params_i__, std::vector<double>& vars__,
                   bool emit_transformed_parameters__ = true,
                   bool emit_generated_quantities__ = true,
                   std::ostream* pstream__ = 0) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "example_model_namespace::write_array";
(void) function__;  // suppress unused var warning

    (void) function__;  // suppress unused var warning

    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    
    try {
      double cont;
      
      current_statement__ = 1;
      cont = in__.scalar();
      current_statement__ = 1;
      cont = stan::math::lub_constrain(cont, -1000, 1000);
      double slope;
      
      current_statement__ = 2;
      slope = in__.scalar();
      current_statement__ = 2;
      slope = stan::math::lub_constrain(slope, -100, 100);
      double amp;
      
      current_statement__ = 3;
      amp = in__.scalar();
      current_statement__ = 3;
      amp = stan::math::lub_constrain(amp, 0, 1000);
      double center;
      
      current_statement__ = 4;
      center = in__.scalar();
      current_statement__ = 4;
      center = stan::math::lub_constrain(center, 0, 10);
      double width;
      
      current_statement__ = 5;
      width = in__.scalar();
      current_statement__ = 5;
      width = stan::math::lub_constrain(width, 0, 10);
      vars__.push_back(cont);
      vars__.push_back(slope);
      vars__.push_back(amp);
      vars__.push_back(center);
      vars__.push_back(width);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array() 
    
  void transform_inits(const stan::io::var_context& context__,
                       std::vector<int>& params_i__,
                       std::vector<double>& vars__, std::ostream* pstream__) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      
      pos__ = 1;
      double cont;
      
      current_statement__ = 1;
      cont = context__.vals_r("cont")[(1 - 1)];
      current_statement__ = 1;
      cont = stan::math::lub_free(cont, -1000, 1000);
      double slope;
      
      current_statement__ = 2;
      slope = context__.vals_r("slope")[(1 - 1)];
      current_statement__ = 2;
      slope = stan::math::lub_free(slope, -100, 100);
      double amp;
      
      current_statement__ = 3;
      amp = context__.vals_r("amp")[(1 - 1)];
      current_statement__ = 3;
      amp = stan::math::lub_free(amp, 0, 1000);
      double center;
      
      current_statement__ = 4;
      center = context__.vals_r("center")[(1 - 1)];
      current_statement__ = 4;
      center = stan::math::lub_free(center, 0, 10);
      double width;
      
      current_statement__ = 5;
      width = context__.vals_r("width")[(1 - 1)];
      current_statement__ = 5;
      width = stan::math::lub_free(width, 0, 10);
      vars__.push_back(cont);
      vars__.push_back(slope);
      vars__.push_back(amp);
      vars__.push_back(center);
      vars__.push_back(width);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits() 
    
  void get_param_names(std::vector<std::string>& names__) const {
    
    names__.resize(0);
    names__.push_back("cont");
    names__.push_back("slope");
    names__.push_back("amp");
    names__.push_back("center");
    names__.push_back("width");
    } // get_param_names() 
    
  void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.resize(0);
    std::vector<size_t> dims__;
    dimss__.push_back(dims__);
    dims__.resize(0);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dimss__.push_back(dims__);
    dims__.resize(0);
    
    } // get_dims() 
    
  void constrained_param_names(std::vector<std::string>& param_names__,
                               bool emit_transformed_parameters__ = true,
                               bool emit_generated_quantities__ = true) const {
    
    param_names__.push_back(std::string() + "cont");
    param_names__.push_back(std::string() + "slope");
    param_names__.push_back(std::string() + "amp");
    param_names__.push_back(std::string() + "center");
    param_names__.push_back(std::string() + "width");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  void unconstrained_param_names(std::vector<std::string>& param_names__,
                                 bool emit_transformed_parameters__ = true,
                                 bool emit_generated_quantities__ = true) const {
    
    param_names__.push_back(std::string() + "cont");
    param_names__.push_back(std::string() + "slope");
    param_names__.push_back(std::string() + "amp");
    param_names__.push_back(std::string() + "center");
    param_names__.push_back(std::string() + "width");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"cont\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"slope\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"amp\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"center\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"width\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"cont\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"slope\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"amp\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"center\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"width\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool emit_transformed_parameters__ = true,
                     bool emit_generated_quantities__ = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng__, params_r_vec, params_i_vec, vars_vec,
          emit_transformed_parameters__, emit_generated_quantities__, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    template <bool propto__, bool jacobian__, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto__,jacobian__,T_>(vec_params_r, vec_params_i, pstream);
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }

};
}
typedef example_model_namespace::example_model stan_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

#endif


