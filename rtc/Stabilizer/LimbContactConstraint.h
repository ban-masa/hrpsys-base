#ifndef LIMBCONTACTCONSTRAINT_H
#define LIMBCONTACTCONSTRAINT_H
#define CONE_DIM 4

enum ConstraintType {
  SURFACE,
  ROPE
};

class LimbContactConstraint
{
public:
  hrp::dmatrix ee_contact_matrix;
  hrp::dmatrix ee_zmp_matrix;
  std::string limb_name;
  //shape info
  std::vector<Eigen::Vector2d> local_foot_vertices;
  std::vector<hrp::Vector3> world_foot_vertices;
  size_t vertices_num;
  //friction info
  Eigen::Vector2d static_friction_coefficients;
  hrp::dmatrix friction_cone_matrix;

  ConstraintType constraint_type;

  LimbContactConstraint (): limb_name("") {}
  LimbContactConstraint (const std::string& _limb_name)
  {
    limb_name = _limb_name;
    Eigen::Vector2d sfc;
    sfc(0) = 0.5;
    sfc(1) = 0.5;
    set_static_friction_coefficients(sfc);
  }
  ~LimbContactConstraint () {}

  //setter
  void set_limb_name (const std::string& _limb_name) { limb_name = _limb_name; };
  void set_constraint_type (ConstraintType _type) { constraint_type = _type; };
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
  void get_limb_name (std::string& _limb_name ) { _limb_name = limb_name; };
  void get_contact_matrix (hrp::dmatrix& mat) { mat = ee_contact_matrix; };
  void get_zmp_matrix (hrp::dmatrix& mat) { mat = ee_zmp_matrix; };
  int get_cone_dim(void) { return CONE_DIM; };
  int get_vertices_num(void) { return vertices_num; };
  int get_state_dim (void)
  {
    int state_dim = 0;
    switch (constraint_type) {
      case SURFACE:
        state_dim = CONE_DIM * vertices_num;
        break;
      case ROPE:
        state_dim = 1;
        break;
      default:
        break;
    }
    return state_dim;
  };
  void get_static_friction_coefficients (Eigen::Vector2d& sfc)
  {
    sfc = static_friction_coefficients;
  };

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

  void calcSurfaceContactMatrix (const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot)
  {
    ee_contact_matrix = hrp::dmatrix::Zero(6, vertices_num * CONE_DIM);
    for (size_t i = 0; i < vertices_num; i++) {
      hrp::dmatrix contact_mat(6, CONE_DIM);
      calcWorldFrictionConeMatrix(contact_mat, i, ee_pos, ee_rot);
      ee_contact_matrix.block(0, i * CONE_DIM, 6, CONE_DIM) = contact_mat;
    }
  }

  void calcSurfaceZMPMatrix (const hrp::Vector3& ee_pos, const hrp::Matrix33& ee_rot, const hrp::Vector3& ref_zmp)
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

  void calcRopeContactMatrix (const hrp::Vector3& ee_pos, const hrp::Vector3& rope_dir)
  {
    ee_contact_matrix = hrp::dmatrix::Zero(6, 1);
    for (size_t i = 0; i < 3; i++) ee_contact_matrix(i, 0) = rope_dir(i);
    hrp::Matrix33 vertex_cross_product;
    vertex_cross_product << 0, -ee_pos(2), ee_pos(1),
                            ee_pos(2), 0, -ee_pos(0),
                            -ee_pos(1), ee_pos(0), 0;
    ee_contact_matrix.block(3, 0, 3, 1) = vertex_cross_product * rope_dir;
  }

  void calcRopeZMPMatrix (void)
  {
    ee_zmp_matrix = hrp::dmatrix::Zero(3, 1);
  }
};
#endif
