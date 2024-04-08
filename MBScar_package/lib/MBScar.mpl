# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#    __  __ ____ ____                                                           #
#   |  \/  | __ ) ___|  ___ __ _ _ __                                           #
#   | |\/| |  _ \___ \ / __/ _` | '__|                                          #
#   | |  | | |_) |__) | (_| (_| | |                                             #
#   |_|  |_|____/____/ \___\__,_|_|                                             #
#                                                                               #
#                                                                               #
#        An integration package for MBSymba                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Matteo Larcher (University of Trento)
#

(*
TODO:
  - add "convenient" reference frame to simplify formulation of symmetrical bodies (see end of page 238 book)
  - add and test constraints handling
*)

MBScar := module()

  description "A package for the generation of the equations of motion of a "
    "multibody system in minimal coordinate and quasi-velocities form. The theory "
    " behind the package is based on the definitions found in the book 'Advanced "
    "Dynamics' by Donald T. Greenwood.\n"
    "This moule extend the functionalities of the MBSymba_r6 module by "
    "R.Lot (c)2003-2011 and M.Massaro (c).";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  # load required packages
  use MBSymba_r6 in

  # local variables
  local m_WarningMode  := true; # default warning mode
  local m_TimeLimit    := 10.0; # default time limit for simplification
  local m_ParallelMode := true; # default parallel mode
  #local m_CacheSize    := 20;   # default cache size for the chache tables

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info := proc()

    description "Print module information.";

    printf(
      "+--------------------------------------------------------------------+\n"
      "| 'MBScar' module information:                                       |\n"
      "| Current version authors:                                           |\n"
      "|   Matteo Larcher.                                                  |\n"
      "+--------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad := proc()

    description "Module load procedure.";

    local lib_base_path, i, types;

    lib_base_path := NULL;
    for i in [libname] do
      if evalb(StringTools:-Search("MBScar", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if evalb(lib_base_path = NULL) then
      error("cannot find 'MBScar' library in the toolbox folder.");
    end if;

    # Define module protected types
    TypeTools:-AddType('MBSymba_BODY',   IsBODY);   protect('MBSymba_BODY');
    TypeTools:-AddType('MBSymba_FORCE',  IsFORCE);  protect('MBSymba_FORCE');
    TypeTools:-AddType('MBSymba_TORQUE', IsTORQUE); protect('MBSymba_TORQUE');
    TypeTools:-AddType('MBSymba_POINT',  IsPOINT);  protect('MBSymba_POINT');
    TypeTools:-AddType('MBSymba_VECTOR', IsVECTOR); protect('MBSymba_VECTOR');

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload := proc()

    description "Module unload procedure.";

    # unprotect module types

    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetModuleOptions := proc(
    {
    WarningMode::{boolean, nothing}   := NULL,
    TimeLimit::{nonnegative, nothing} := NULL
    }, $)

    description "Set the module options: warning mode <WarningMode>, time limit "
      "<TimeLimit>.";

    if evalb(WarningMode <> NULL) then
      m_WarningMode := WarningMode;
    end if;

    if evalb(TimeLimit <> NULL) then
      m_TimeLimit := TimeLimit;
    end if;

    return NULL;
  end proc: # SetModuleOptions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsBODY := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of BODY type.";

  return type(var, table) and evalb(var[parse("obj")] = BODY);
  end proc: # IsBODY

  # - - - - - - - - - - - - - - - -

  export IsFORCE := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of FORCE type.";

  return type(var, table) and evalb(var[parse("obj")] = FORCE);
  end proc: # IsFORCE

  # - - - - - - - - - - - - - - - -

  export IsTORQUE := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of TORQUE type.";

  return type(var, table) and evalb(var[parse("obj")] = TORQUE);
  end proc: # IsTORQUE

  # - - - - - - - - - - - - - - - -

  export IsPOINT := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of POINT type.";

  return type(var, table) and evalb(var[parse("obj")] = POINT);
  end proc: # IsPOINT

  # - - - - - - - - - - - - - - - -

  export IsVECTOR := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of VECTOR type.";

  return type(var, table) and evalb(var[parse("obj")] = VECTOR);
  end proc: # IsVECTOR

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsEqual := proc(
  obj1::anything,
  obj2::anything,
  $)::boolean;

  description "Check if the two objects <obj1> and <obj2> are equal.";
  option remember;
  local fields;

  # check if the two objects are tables (return false if not)
  if not(type(obj1, table) and type(obj2, table)) then
    return false;
  end if;

  fields := lhs~(op(op(obj1)));

  try
    map(x -> ArrayTools:-IsEqual(obj1[x], obj2[x]), fields);
    return not(has(%, false));
  catch:
    return false;
  end try;
  end proc: # IsEqual



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Simplify := proc(
  var::anything,
  opt::anything := NULL,
  $)::anything;

    description "Try to simplify an algebraic expression <var> with optional "
      "simplification options <opt>. The simplification is performed within "
      "the internal or indexed time limit.";
    option remember;

    try
      return timelimit(
        `if`(procname::indexed, op(procname), m_TimeLimit),
        (simplify(var, opt) assuming real)
      );
    catch "time expired":
      if m_WarningMode then
        WARNING("exceeded time limit, raw solutions is returned.");
      end if;
      return var;
    catch "division by zero":
      error("division by zero detected.");
      return var;
    catch:
      error("something went wrong, last exception '%1'.", lastexception);
      return var;
    end try:
  end proc: # Simplify

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Spy := proc(
  A::Matrix,
  {
    row_labels::{list, nothing} := [seq(i, i=1..LinearAlgebra:-RowDimension(A))],
    col_labels::{list, nothing} := [seq(i, i=1..LinearAlgebra:-ColumnDimension(A))]
  },
  $)::anything;

    description "Plot of non-zero values of matrix <A>. The optional arguments "
      "<row_labels> and <col_labels> are used to label the rows and columns of "
      "the matrix.";
    local opt, i;

    opt := {
      axis[1]=[tickmarks=[[seq(i, i=1..LinearAlgebra:-ColumnDimension(A))] =~ col_labels, rotation = Pi/3], gridlines=[[seq(i+0.5, i = 1..LinearAlgebra:-ColumnDimension(A)-1)]]],
      axis[2]=[tickmarks=[-[seq(i, i=1..LinearAlgebra:-RowDimension(A))] =~ row_labels], gridlines=[[seq(-i-0.5, i = 1..LinearAlgebra:-RowDimension(A)-1)]]]
      };

    return plots:-sparsematrixplot(A, 'matrixview', op(opt));
  end proc: # Spy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export kinetic_energy_qv := proc(
  bodies::{list(MBSymba_BODY),set(MBSymba_BODY)},
  q_vars::list(scalar),
  qv_eqns::list({algebraic, `=`}),
  $)::algebraic;

  description "Compute the kinetic energy of the system in terms of the quasi-velocities.";
  option remember;
  local body, v_vec, omega_vec, VV, MM, K;

  # if m_WarningMode then
  #   if nops(bodies) > m_CacheSize then
  #     WARNING("number of bodies exceeds tne chache size, consider setting a larger "
  #     "<m_CacheSize> in module options to increase performances");
  #   end if;
  # end if;

  K := 0;

  userinfo(3, kinetic_energy_qv, "computing kinetic energy");
  for body in bodies do
    # linear velocity vector v_vec
    v_vec := linear_velocity_qv(body, q_vars, qv_eqns);
    # angular velocity vector omega_vec
    omega_vec := angular_velocity_qv(body, q_vars, qv_eqns);

    # velocities vector VV
    VV := [op(v_vec), op(omega_vec)];
    userinfo(5, kinetic_energy_qv, "velocity vector", print(VV));

    # mass matrix
    MM := LinearAlgebra:-DiagonalMatrix([Matrix(3,3,shape=diagonal,fill=body[parse("mass")]), body[parse("inertia")]]);
    userinfo(5, kinetic_energy_qv, "mass matrix", print(MM));

    # compute kinetic energy according to Koenig’s theorem and add to the total kinetic energy
    K := K + 1/2 * (<VV>^%T . MM . <VV>)
  end do;
  K := Simplify(K);
  userinfo(5, kinetic_energy_qv, "kinetic energy", print(K));
  userinfo(3, kinetic_energy_qv, "computing kinetic energy -- DONE");

  return K;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export linear_velocity_qv := proc(
  body::MBSymba_BODY,
  q_vars::list(scalar),
  qv_eqns::list({algebraic, `=`}),
  $)

  description "Compute the linear velocity of the body terms of the quasi-velocities.";
  option remember;
  local v_vec;

  userinfo(3, linear_velocity_qv, "computing linear velocity");
  v_vec := simplify(simplify([MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-velocity(MBSymba_r6_kinematics:-origin(body[parse("frame")])),body[parse("frame")])], qv_eqns, diff(q_vars,t)),trig);
  if has(v_vec, diff(q_vars,t)) then
    error "linear velocity contains time derivatives of the generalized coordinates";
  end if;
  userinfo(5, linear_velocity_qv, "linear velocity vector", print(v_vec));
  userinfo(3, linear_velocity_qv, "computing linear velocity -- DONE");

  return Simplify(v_vec);
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export angular_velocity_qv := proc(
  body::MBSymba_BODY,
  q_vars::list(scalar),
  qv_eqns::list({algebraic, `=`}),
  $)

  description "Compute the angular velocity of the body terms of the quasi-velocities.";
  option remember;
  local omega_vec;

  userinfo(3, angular_velocity_qv, "computing angular velocity");
  omega_vec := simplify(simplify([MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-angular_velocity(body[parse("frame")]),body[parse("frame")])], qv_eqns, diff(q_vars,t)),trig);
  if has(omega_vec, diff(q_vars,t)) then
    error "angular velocity contains time derivatives of the generalized coordinates";
  end if;
  userinfo(5, angular_velocity_qv, "angular velocity vector", print(omega_vec));
  userinfo(3, angular_velocity_qv, "computing angular velocity -- DONE");

  return Simplify(omega_vec);
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export linear_momentum_qv := proc(
  body::MBSymba_BODY,
  q_vars::list(scalar),
  qv_eqns::list({algebraic, `=`}),
  $)::algebraic;

  description "Compute the linear momentum of the body terms of the quasi-velocities.";
  option remember;
  local p, v_vec;

  userinfo(3, linear_momentum_qv, "computing linear momentum");

  # linear velocity vector v_vec
  v_vec := linear_velocity_qv(body, q_vars, qv_eqns);

  # linear momentum
  p := body[parse("mass")] *~ v_vec;
  userinfo(5, linear_momentum_qv, "linear momentum", print(p));
  userinfo(3, linear_momentum_qv, "computing linear momentum -- DONE");

  return Simplify(p);
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export angular_momentum_qv := proc(
  body::MBSymba_BODY,
  q_vars::list(scalar),
  qv_eqns::list({algebraic, `=`}),
  $)::algebraic;

  description "Compute the angular momentum of the body terms of the quasi-velocities.";
  option remember;
  local H, omega_vec;

  userinfo(3, angular_momentum_qv, "computing angular momentum");

  # angular velocity vector omega_vec
  omega_vec := angular_velocity_qv(body, q_vars, qv_eqns);

  # angular momentum
  H := convert(body[parse("inertia")] . <omega_vec>, list);
  userinfo(5, angular_momentum_qv, "angular momentum", print(H));
  userinfo(3, angular_momentum_qv, "computing angular momentum -- DONE");

  return Simplify(H);
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export fundamental_equations := proc(
  mbo::{list({MBSymba_BODY,MBSymba_FORCE,MBSymba_TORQUE}),set({MBSymba_BODY,MBSymba_FORCE,MBSymba_TORQUE})},
  q_vars::list(scalar),
  u_vars::list(scalar),
  qv_eqns::list({algebraic, `=`}),
  $)::list(algebraic);

  description "Compute the fundamental equations of motion of the system.";
  option remember;
  local bodies, forces, obj, eqns, body, force, T, p, H, gamma, gamma_dot, beta, beta_dot,
    v, omega, u, Q, i, j;

  # if m_WarningMode then
  #   if nops(bodies) > m_CacheSize then
  #     WARNING("number of bodies exceeds tne chache size, consider setting a larger "
  #     "<m_CacheSize> in module options to increase performances");
  #   end if;
  # end if;

  userinfo(3, fundamental_equations, "computing fundamental equations of motion");

  # separate bodies and forces
  bodies := [];
  forces := [];
  for obj in mbo do
    if obj[parse("obj")] = BODY then
      bodies := [op(bodies), copy(obj, deep)]; # copy the object to avoid modifications to the original object
    elif obj[parse("obj")] = FORCE or obj[parse("obj")] = TORQUE then
      forces := [op(forces), copy(obj, deep)]; # copy the object to avoid modifications to the original object
    end if;
  end do;

  # create bodies gravity forces
  userinfo(3, fundamental_equations, "creating gravity forces");
  for body in bodies do
    forces := [op(forces), MBSymba_r6_dynamics:-make_FORCE(MBSymba_r6_kinematics:-make_VECTOR(_gravity[parse("frame")], op([MBSymba_r6_kinematics:-comp_XYZ(_gravity)] *~ body[parse("mass")])), MBSymba_r6_kinematics:-origin(body[parse("frame")]), body)];
    userinfo(5, fundamental_equations, "gravity force", print(show(forces[-1])));
  end do;
  userinfo(3, fundamental_equations, "creating gravity forces -- DONE");

  # project forces on bodies origin (origin reference point can be different from the body center of mass)
  userinfo(3, fundamental_equations, "projecting forces in bodies center of mass");
  for force in forces do
    if IsFORCE(force) then
      # add the additional torque to the forces list (FIXME: MBSymba_r6_kinematics:-cross_product have performance issues, consider to implement a custom cross/dot product)
      MBSymba_r6_kinematics:-cross_prod(MBSymba_r6_kinematics:-join_points(MBSymba_r6_kinematics:-origin(force[parse("acting")][parse("frame")]), force[parse("applied")]), force);
      forces := [op(forces), %];
      userinfo(5, fundamental_equations, "translational torque", print(show(forces[-1])));
      # create "acting" torque on reacting body
      if has(lhs~(op(op(force))), parse("reacting")) then
        MBSymba_r6_kinematics:-cross_prod(MBSymba_r6_kinematics:-join_points(MBSymba_r6_kinematics:-origin(force[parse("reacting")][parse("frame")]), force[parse("applied")]), force);
        %[parse("comps")] := -%[parse("comps")];
        forces := [op(forces), %];
        userinfo(5, fundamental_equations, "translational torque (reacting body)", print(show(forces[-1])));
        force[parse("reacting")] := NULL;
      end if;
      # change application point to body origin (reference point)
      force[parse("applied")] := evala(MBSymba_r6_kinematics:-origin(force[parse("acting")][parse("frame")]));
      userinfo(5, fundamental_equations, "projected force", print(show(force)));
      # create "acting" force on reacting body
      if has(lhs~(op(op(force))), parse("reacting")) then
        MBSymba_r6_dynamics:-make_FORCE(MBSymba_r6_kinematics:-make_VECTOR(-force[parse("comps")], force[parse("frame")]), MBSymba_r6_kinematics:-origin(force[parse("reacting")][parse("frame")]), force[parse("reacting")]);
        forces := [op(forces), %];
        userinfo(5, fundamental_equations, "projected force (reacting body)", print(show(forces[-1])));
      end if;
    end if;
  end do;
  userinfo(3, fundamental_equations, "projecting forces in bodies center of mass -- DONE");

  # compute bodies velocities vectors
  userinfo(3, fundamental_equations, "computing bodies velocities vectors");
  v := [seq(0, i=1..nops(bodies))]; i := 1;
  for body in bodies do
    v[i] := linear_velocity_qv(body, q_vars, qv_eqns);
    userinfo(4, fundamental_equations, "linear velocity", i, " = ", print(v[i]));
    i := i + 1;
  end do;
  userinfo(3, fundamental_equations, "computing bodies velocities vectors -- DONE");

  # compute bodies angular velocity vectors
  userinfo(3, fundamental_equations, "computing bodies angular velocity vectors");
  omega := [seq(0, i=1..nops(bodies))]; i := 1;
  for body in bodies do
    omega[i] := angular_velocity_qv(body, q_vars, qv_eqns);
    userinfo(4, fundamental_equations, "angular velocity", i, " = ", print(omega[i]));
    i := i + 1;
  end do;
  userinfo(3, fundamental_equations, "computing bodies angular velocity vectors -- DONE");

  # compute gammas and their derivatives
  userinfo(3, fundamental_equations, "computing gammas");
  gamma := [seq([seq(0, j=1..nops(u_vars))], i=1..nops(bodies))];
  gamma_dot := [seq([seq(0, j=1..nops(u_vars))], i=1..nops(bodies))];
  i := 1; j := 1;
  for body in bodies do
    for u in u_vars do
      gamma[i,j] := Physics:-diff(v[i], u); #(4.223 book)
      userinfo(4, fundamental_equations, "gamma", i, j, " = ", print(gamma[i,j]));
      gamma_dot[i,j] := diff(gamma[i,j],t) +~ convert(LinearAlgebra:-CrossProduct(omega[i], gamma[i,j]), list); # Poisson's formula
      userinfo(4, fundamental_equations, "gamma_dot", i, j, " = ", print(gamma_dot[i,j]));
      j := j + 1;
    end do;
    j := 1;
    i := i + 1;
  end do;
  userinfo(3, fundamental_equations, "computing gammas -- DONE");

  # compute betas and their derivatives
  userinfo(3, fundamental_equations, "computing betas");
  beta := [seq([seq(0, j=1..nops(u_vars))], i=1..nops(bodies))];
  beta_dot := [seq([seq(0, j=1..nops(u_vars))], i=1..nops(bodies))];
  i := 1; j := 1;
  for body in bodies do
    for u in u_vars do
      beta[i,j] := Physics:-diff(omega[i], u); #(4.223 book)
      userinfo(4, fundamental_equations, "beta", i, j, " = ", print(beta[i,j]));
      beta_dot[i,j] := diff(beta[i,j],t) +~ convert(LinearAlgebra:-CrossProduct(omega[i], beta[i,j]), list); # Poisson's formula
      userinfo(4, fundamental_equations, "beta_dot", i, j, " = ", print(beta_dot[i,j]));
      j := j + 1;
    end do;
    j := 1;
    i := i + 1;
  end do;
  userinfo(3, fundamental_equations, "computing betas -- DONE");

  # compute generalized forces
  userinfo(3, fundamental_equations, "computing generalized forces");
  Q := [seq(0, j=1..nops(u_vars))]; i := 1; j := 1;
  for u in u_vars do
    for body in bodies do
        for force in forces do
          if IsFORCE(force) then
            if IsEqual(force[parse("acting")], body) then
              Q[j] := Q[j] + <MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-project(force, body[parse("frame")]))>^%T . <gamma[i,j]>;
            elif IsEqual(force[parse("reacting")], body) then # NOTE: no force should have a reacting body after projection
              Q[j] := Q[j] + <-MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-project(force, body[parse("frame")]))>^%T . <gamma[i,j]>;
            end if;
          elif IsTORQUE(force) then
            if IsEqual(force[parse("acting")], body) then
              Q[j] := Q[j] + <MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-project(force, body[parse("frame")]))>^%T . <beta[i,j]>;
            elif IsEqual(force[parse("reacting")], body) then
              Q[j] := Q[j] + <-MBSymba_r6_kinematics:-comp_XYZ(MBSymba_r6_kinematics:-project(force, body[parse("frame")]))>^%T . <beta[i,j]>;
            end if;
          end if;
          userinfo(4, fundamental_equations, "generalized force", j, " = ", print(Q[j]));
        end do;
        i := i + 1;
    end do;
    i := 1;
    j := j + 1;
  end do;
  userinfo(3, fundamental_equations, "computing generalized forces -- DONE");

  # compute bodies linear momentum
  userinfo(3, fundamental_equations, "computing bodies linear momentum");
  p := [seq(0, i=1..nops(bodies))]; i := 1;
  for body in bodies do
    p[i] := linear_momentum_qv(body, q_vars, qv_eqns);
    userinfo(4, fundamental_equations, "linear momentum", i, " = ", print(p[i]));
    i := i + 1;
  end do;
  userinfo(3, fundamental_equations, "computing bodies linear momentum -- DONE");

  # compute bodies linear momentum
  userinfo(3, fundamental_equations, "computing bodies angular momentum");
  H := [seq(0, i=1..nops(bodies))]; i := 1;
  for body in bodies do
    H[i] := angular_momentum_qv(body, q_vars, qv_eqns);
    userinfo(4, fundamental_equations, "angular momentum", i, " = ", print(H[i]));
    i := i + 1;
  end do;
  userinfo(3, fundamental_equations, "computing bodies angular momentum -- DONE");

  # compute kinetic energy
  userinfo(3, fundamental_equations, "computing kinetic energy");
  T := kinetic_energy_qv(bodies, q_vars, QV_eq);
  userinfo(4, fundamental_equations, "kinetic energy", print(T));
  userinfo(3, fundamental_equations, "computing kinetic energy -- DONE");

  # compute fundamental equation of motion
  userinfo(3, fundamental_equations, "computing fundamental equations of motion");
  eqns := [seq(0, j=1..nops(u_vars))]; j := 1; i := 1;
  for u in u_vars do
    simplify(simplify(diff(Physics:-diff(T,u),t), qv_eqns, diff(q_vars,t)),trig); # d(d(T)/du)/dt
    if has(%, diff(q_vars,t)) then
      error "linear velocity contains time derivatives of the generalized coordinates";
    end if;
    eqns[j] := % - Q[j];
    for body in bodies do
      # fundamental equation of motion (4.227 book)
      eqns[j] := eqns[j] - <p[i]>^%T . <gamma_dot[i,j]> - <H[i]>^%T . <beta_dot[i,j]>;
      i := i + 1;
    end do;
    userinfo(4, fundamental_equations, "equation of motion", j, " = ", print(eqns[j]));
    i := 1;
    j := j + 1;
  end do;

  userinfo(3, fundamental_equations, "computing fundamental equations of motion -- DONE");

  return Simplify(eqns);
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end use:
end module: # MBscar_EQM

# That's all folks!