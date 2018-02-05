void set_matrix_with_rope (
    hrp::dmatrix& Phimat, const std::vector<bool>& is_contact_list,
    const hrp::Vector3& rope_dir, const hrp::Vector3& hand_zmp, const size_t grasping_hand_num)
{
  size_t state_dim = 0;
  size_t ee_contact_num = 0;
  for (size_t i = 0; i < is_contact_list.size(); i++) {
    if (is_contact_list[i]) {
      state_dim += ee_matrix_list[i].cols();
      ee_contact_num++;
    }
  }
  if (grasping_hand_num > 0) state_dim++; //For rope state variable
  Phimat = hrp::dmatrix::Zero(6 + state_dim, state_dim);
  size_t col_idx = 0;
  for (size_t i = 0; i < is_contact_list.size(); i++) {
    if (is_contact_list[i]) {
      Phimat.block(0, col_idx, 6, ee_matrix_list[i].cols()) = ee_matrix_list[i];
      col_idx += ee_matrix_list[i].cols();
    }
  }

  if (grasping_hand_num > 0) {
    for (int i = 0; i < 3; i++) Phimat(i, state_dim - 1) = rope_dir(i);
    hrp::Vector3 temp_vec = hand_zmp.cross(rope_dir);
    for (int i = 0; i < 3; i++) Phimat(i + 3, state_dim - 1) = temp_vec(i);
  }
  Phimat.block(6, 0, state_dim, state_dim) = hrp::dmatrix::Identity(state_dim, state_dim);
}

void calc_hand_zmp_and_rope_dir(hrp::Vector3& hand_zmp, hrp::Vector3& rope_dir,
                                const std::vector<hrp::Vector3>& act_hand_force,
                                const std::vector<hrp::Vector3>& act_hand_moment,
                                const std::vector<hrp::Vector3>& hands_pos,
                                const std::vector<hrp::Matrix33>& hands_rot,
                                const std::vector<bool>& hand_contact_list,
                                const size_t grasping_hand_num
                                )
{
  hrp::Vector3 total_force = hrp::Vector3::Zero();
  hrp::Vector3 total_moment = hrp::Vector3::Zero();
  for (size_t i = 0; i < hand_contact_list.size(); i++) {
    if (hand_contact_list[i]) {
      total_force = total_force + act_hand_force[i];
      total_moment = total_moment + hands_pos[i].cross(act_hand_force[i]) + act_hand_moment[i];
    }
  }
  if (grasping_hand_num == 0) {
    return;
  } else if (grasping_hand_num == 1) {
    for (size_t i = 0; i < hand_contact_list.size(); i++) {
      if (hand_contact_list[i]) hand_zmp = hands_pos[i];
    }
  } else {
    hrp::Matrix33 temp_mat;
    temp_mat << 0, total_force(2), -total_force(1),
             -total_force(2), 0, total_force(0),
             total_force(1), -total_force(0), 0;
    hand_zmp = temp_mat.inverse() * total_moment;
  }
  rope_dir = total_force.normalized();
}

void calc_hand_ref_wrench(const double rope_tension,
                          std::vector<hrp::Vector3>& ref_hand_force,
                          std::vector<hrp::Vector3>& ref_hand_moment,
                          const hrp::Vector3& hand_zmp,
                          const hrp::Vector3& rope_dir,
                          const std::vector<hrp::Vector3>& hands_pos,
                          const std::vector<hrp::Matrix33>& hands_rot,
                          const std::vector<bool>& hand_contact_list,
                          const size_t grasping_hand_num
    )
{
  hrp::Vector3 rope_force = rope_dir * rope_tension;
  ref_hand_moment[0] = hrp::Vector3::Zero();
  ref_hand_moment[1] = hrp::Vector3::Zero();
  if (grasping_hand_num == 0) {
    ref_hand_force[0] = hrp::Vector3::Zero();
    ref_hand_force[1] = hrp::Vector3::Zero();
    return;
  } else if (grasping_hand_num == 1) {
    for (size_t i = 0; i < hand_contact_list.size(); i++) {
      if (hand_contact_list[i]) {
        ref_hand_force[i] = rope_force;
      } else {
        ref_hand_force[i] = hrp::Vector3::Zero();
      }
    }
  } else if (grasping_hand_num == 2) {
    hrp::Matrix33 cross_mat1, cross_mat2;
    cross_mat1 << 0, -hands_pos[0](2), hands_pos[0](1),
                  hands_pos[0](2), 0, -hands_pos[0](0),
                  -hands_pos[0](1), hands_pos[0](0), 0;
    cross_mat2 << 0, -hands_pos[1](2), hands_pos[1](1),
                  hands_pos[1](2), 0, -hands_pos[1](0),
                  -hands_pos[1](1), hands_pos[1](0), 0;
    ref_hand_force[0] = -(cross_mat1 - cross_mat2).inverse() * cross_mat2 * rope_force;
    ref_hand_force[1] = rope_force - ref_hand_force[0];
  }
}

void distributeZMPToForceMomentsQPAlllimbs (
    std::vector<hrp::Vector3>& ref_foot_force, std::vector<hrp::Vector3>& ref_foot_moment,
    std::vector<hrp::Vector3>& ref_hand_force, std::vector<hrp::Vector3>& ref_hand_moment,
    const std::vector<hrp::Vector3>& act_hand_force, const std::vector<hrp::Vector3>& act_hand_moment,
    const std::vector<hrp::Vector3>& ee_pos,
    const std::vector<hrp::Vector3>& hands_pos,
    const std::vector<hrp::Matrix33>& ee_rot,
    const std::vector<hrp::Matrix33>& hands_rot,
    const std::vector<std::string>& ee_name,
    const std::vector<std::string>& hands_name,
    const hrp::Vector3& ref_zmp, const hrp::Vector3& ref_cog,
    const double total_mass, const double gravity,
    const double total_fz,
    const std::vector<bool>& is_contact_list, const std::vector<bool>& hand_contact_list,
    const bool printp = true, const std::string& print_str = "")
{
  size_t ee_num = ee_name.size();
  size_t contact_ee_num = 0;
  for (size_t i = 0; i < is_contact_list.size(); i++) {//足の着地状態をカウント
    if (is_contact_list[i]) contact_ee_num++;
  }
  if (contact_ee_num == 0) {//両足とも接地していないならreturn
    for (size_t i = 0; i < ee_num; i++) {
      ref_foot_force[i] = hrp::Vector3(0, 0, total_fz);
      ref_foot_moment[i] = hrp::Vector3::Zero();
    }
    return;
  }

  size_t grasping_hand_num = 0;//把持しているハンドの数
  for (size_t i = 0; i < hand_contact_list.size(); i++) {
    if (hand_contact_list[i]) {
      grasping_hand_num++;
    }
  }
  
  hrp::dvector total_wrench(6);//目標wrench
  double omega2 = total_fz / (ref_cog(2) - ref_zmp(2));
  total_wrench.head(3) = hrp::Vector3(omega2 * (ref_cog(0) - ref_zmp(0)), omega2 * (ref_cog(1) - ref_zmp(1)), total_fz);
  total_wrench.tail(3) = ref_cog.cross((hrp::Vector3)total_wrench.head(3));

  size_t cone_dim = 4;//摩擦錘の底面の頂点数
  size_t state_dim = 0;//状態変数のサイズをカウント
  for (size_t i = 0; i < ee_num; i++) {
    if (is_contact_list[i]) {
      state_dim += cone_dim * fs.get_vertices_num(i);
    }
  }
  if (grasping_hand_num > 0) state_dim++;//把持ハンドがあるならば+1

  hrp::Vector3 hand_zmp;
  hrp::Vector3 rope_dir;
  calc_hand_zmp_and_rope_dir(hand_zmp, rope_dir,
                             act_hand_force, act_hand_moment,
                             hands_pos, hands_rot, hand_contact_list, grasping_hand_num);
  hrp::dmatrix Phimat = hrp::dmatrix::Zero(6 + state_dim, state_dim);
  hrp::dvector xivec = hrp::dvector::Zero(6 + state_dim);
  hrp::dmatrix Wmat = hrp::dmatrix::Zero(6 + state_dim, 6 + state_dim);
  hrp::dmatrix Hmat = hrp::dmatrix::Zero(state_dim, state_dim);
  hrp::dvector gvec = hrp::dvector::Zero(state_dim);

  calcAllContactMatrix(ee_pos, ee_rot, ee_name, cone_dim);
  set_matrix_with_rope(Phimat, is_contact_list, rope_dir, hand_zmp, grasping_hand_num);
  //TODO: make mode to adapt to force, moment weight vector

  for (size_t i = 0; i < 6; i++) {
    xivec(i) = total_wrench(i);
  }
  //TODO: set correct param to Wmat
  for (size_t i = 0; i < state_dim + 6; i++) {
    if (i < 6) Wmat(i, i) = weight_param_for_qp_weight_matrix[i];
    else Wmat(i, i) = weight_param_for_qp_weight_matrix[6];
  }
  Hmat = Phimat.transpose() * Wmat * Phimat;
  gvec = -xivec.transpose() * Wmat * Phimat;

  std::vector<hrp::dvector> rho_vec;
  for (size_t i = 0; i < ee_num; i++) {
    if (is_contact_list[i]) {
      rho_vec.push_back(hrp::dvector(cone_dim * fs.get_vertices_num(i)));
    }
  }

  solveForceMomentQPOASES(rho_vec, state_dim, contact_ee_num, Hmat, gvec);
  hrp::dvector tmpv(6);
  size_t tmp_idx = 0;
  for (size_t i = 0; i < ee_num; i++) {
    if (is_contact_list[i]) {
      tmpv = ee_matrix_list[i] * rho_vec[tmp_idx];
      tmp_idx++;
      ref_foot_force[i] = hrp::Vector3(tmpv(0), tmpv(1), tmpv(2));
      ref_foot_moment[i] = hrp::Vector3(tmpv(3), tmpv(4), tmpv(5)) - ee_pos[i].cross(ref_foot_force[i]);
    } else {
      ref_foot_force[i] = hrp::Vector3::Zero();
      ref_foot_moment[i] = hrp::Vector3::Zero();
    }
  }
  calc_hand_ref_wrench(rope_tension,
                       ref_hand_force,
                       ref_hand_moment,
                       hand_zmp,
                       rope_dir,
                       hands_pos,
                       hands_rot,
                       hand_contact_list,
                       grasping_hand_num);
                       
}
