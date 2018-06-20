#define CONE_DIM 4

enum LimbType {
  RLEG,
  LLEG,
  RARM,
  LARM
};

class LimbContactConstraintBase
{
protected:
  hrp::dmatrix ee_contact_matrix;
  hrp::dmatrix ee_zmp_matrix;
  std::string limb_name;
public:
  LimbContactConstraintBase (): limb_name("") {}
  LimbContactConstraintBase (const std::string& _limb_name) { limb_name = _limb_name; }
  ~LimbContactConstraintBase () {}

  //setter
  virtual void set_limb_name (const std::string& _limb_name) { limb_name = _limb_name; };
  virtual void set_foot_vertices (const std::vector<Eigen::Vector2d>& fvs) { std::cerr << "set support polygon vertices" << std::endl; };
  virtual void set_static_friction_coefficients (const Eigen::Vector2d& sfc) {};
  //getter
  virtual void get_limb_name (std::string& _limb_name ) { _limb_name = limb_name; };
  virtual void get_contact_matrix (hrp::dmatrix& mat) { mat = ee_contact_matrix; };
  virtual void get_zmp_matrix (hrp::dmatrix& mat) { mat = ee_zmp_matrix; };
  virtual int get_state_dim (void) { return 0; };
  virtual void get_static_friction_coefficients(Eigen::Vector2d& sfc) {};
  //print
  virtual void print_info(void) {};

  virtual void calcContactMatrix(void) { std::cerr << "Implement ContactMatrix method" << std::endl; };
  virtual void calcContactMatrix (const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot) {};
  virtual void calcContactMatrix (const hrp::Vector3& ee_pos, const hrp::Vector3& rope_dir) {};
  virtual void calcZMPMatrix(void) { std::cerr << "Implement ZMPMatrix method" << std::endl; };
  virtual void calcZMPMatrix (const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot, const hrp::Vector3& ref_zmp) {};
};

class SimpleLimbContactConstraint: public LimbContactConstraintBase
{
protected:
  //shape info
  std::vector<Eigen::Vector2d> local_foot_vertices;
  std::vector<hrp::Vector3> world_foot_vertices;
  size_t vertices_num;
  //friction info
  Eigen::Vector2d static_friction_coefficients;
  hrp::dmatrix friction_cone_matrix;
public:
  SimpleLimbContactConstraint (const std::string& _limb_name) : LimbContactConstraintBase(_limb_name), vertices_num(4)
  {
    Eigen::Vector2d sfc;
    sfc(0) = 0.5;
    sfc(1) = 0.5;
    set_static_friction_coefficients(sfc);
  }
  ~SimpleLimbContactConstraint () {}
  //setter
  void set_static_friction_coefficients (const Eigen::Vector2d& sfc)
  {
    static_friction_coefficients = sfc;
    calcFrictionConeMatrix();
  };
  void set_foot_vertices(const std::vector<Eigen::Vector2d>& fvs)
  {
    local_foot_vertices = fvs;
    vertices_num = fvs.size();
    world_foot_vertices.clear();
    for (size_t i = 0; i < vertices_num; i++) {
      world_foot_vertices.push_back(hrp::Vector3::Zero());
    }
  };
  //getter
  virtual void get_static_friction_coefficients (Eigen::Vector2d& sfc) { sfc = static_friction_coefficients; };
  int get_cone_dim(void) { return CONE_DIM; };
  int get_vertices_num(void) { return vertices_num; };
  virtual int get_state_dim (void) { return CONE_DIM * vertices_num; };
  //print
  virtual void print_info (void)
  {
    std::cerr << limb_name << " world vertices" << std::endl;
    for (size_t i = 0; i < world_foot_vertices.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        std::cerr << world_foot_vertices[i](j) << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
  
  //calc friction cone matrix
  void calcFrictionConeMatrix(void)
  {
    friction_cone_matrix = hrp::dmatrix::Zero(3, CONE_DIM);
    for (size_t i = 0; i < CONE_DIM; i++) {
      double theta = i * 2.0 * M_PI / CONE_DIM;
      friction_cone_matrix(0, i) = static_friction_coefficients(0) * std::cos(theta);
      friction_cone_matrix(1, i) = static_friction_coefficients(1) * std::sin(theta);
      friction_cone_matrix(2, i) = 1.0;
    }
  }

  void calcWorldFootVertices (const size_t& idx, const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot)
  {
    Eigen::Vector2d& lfv = local_foot_vertices[idx];
    world_foot_vertices[idx] = ee_pos + ee_rot * hrp::Vector3(lfv(0), lfv(1), 0.0);
  }

  void calcWorldFrictionConeMatrix (hrp::dmatrix& ret_mat, const size_t& idx, const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot)
  {
    calcWorldFootVertices(idx, ee_pos, ee_rot);
    ret_mat = hrp::dmatrix(6, CONE_DIM);
    hrp::dmatrix world_cone_matrix = ee_rot * friction_cone_matrix;
    ret_mat.block(0, 0, 3, CONE_DIM) = world_cone_matrix;
    hrp::Vector3& world_vertex = world_foot_vertices[idx];
    Eigen::Matrix3d vertex_cross_product;
    vertex_cross_product << 0, -world_vertex(2), world_vertex(1),
                            world_vertex(2), 0, -world_vertex(0),
                            -world_vertex(1), world_vertex(0), 0;
    ret_mat.block(3, 0, 3, CONE_DIM) = vertex_cross_product * world_cone_matrix;
  }

  virtual void calcContactMatrix (const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot)
  {
    ee_contact_matrix = hrp::dmatrix::Zero(6, vertices_num * CONE_DIM);
    for (size_t i = 0; i < vertices_num; i++) {
      hrp::dmatrix contact_mat(6, CONE_DIM);
      calcWorldFrictionConeMatrix(contact_mat, i, ee_pos, ee_rot);
      ee_contact_matrix.block(0, i * CONE_DIM, 6, CONE_DIM) = contact_mat;
    }
  }

  virtual void calcZMPMatrix (const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot, const hrp::Vector3& ref_zmp)
  {
    ee_zmp_matrix = hrp::dmatrix::Zero(3, vertices_num * CONE_DIM);
    for (size_t i = 0; i < vertices_num; i++) {
      calcWorldFootVertices(i, ee_pos, ee_rot);
      hrp::Vector3 vertex = world_foot_vertices[i] - ref_zmp;
      hrp::Matrix33 vertex_cross_product;
      vertex_cross_product << 0, -vertex(2), vertex(1),
                              vertex(2), 0, -vertex(0),
                              -vertex(1), vertex(0), 0;
      ee_zmp_matrix.block(0, i * CONE_DIM, 3, CONE_DIM) = vertex_cross_product * ee_rot * friction_cone_matrix;
    }
  }
};

class RopeGraspHandContactConstraint: public LimbContactConstraintBase
{
public:
  RopeGraspHandContactConstraint (const std::string& _limb_name) : LimbContactConstraintBase(_limb_name) {}
  ~RopeGraspHandContactConstraint () {}
  virtual int get_state_dim (void) { return 1; };
  virtual void calcContactMatrix (const hrp::Vector3& ee_pos, const hrp::Vector3& rope_dir)
  {
    ee_contact_matrix = hrp::dmatrix::Zero(6, 1);
    for (size_t i = 0; i < 3; i++) ee_contact_matrix(i, 0) = rope_dir(i);
    hrp::Matrix33 vertex_cross_product;
    vertex_cross_product << 0, -ee_pos(2), ee_pos(1),
                            ee_pos(2), 0, -ee_pos(0),
                            -ee_pos(1), ee_pos(0), 0;
    ee_contact_matrix.block(3, 0, 3, 1) = vertex_cross_product * rope_dir;
  }

  virtual void calcZMPMatrix (void)
  {
    ee_zmp_matrix = hrp::dmatrix::Zero(3, 1);
  }
};
