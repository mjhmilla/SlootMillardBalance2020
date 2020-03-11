return {
  configuration = {
    axis_front = { 1, 0, 0,},
    axis_right = { 0, -1, 0,},
    axis_up = { 0, 0, 1,},
  },
  constraint_sets = {},
  frames = {
    {
      body = {
        com = { 0, 0, 0.084231,},
        inertia = {
          { 0.034067, 0, 0,},
          { 0, 0.029363, 0,},
          { 0, 0, 0.03582,},
        },
        mass = 6.6091,
      },
      joint = {
        { 0, 0, 0, 1, 0, 0,},
        { 0, 0, 0, 0, 1, 0,},
        { 0, 0, 0, 0, 0, 1,},
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, 0,},
      },
      markers = {
        L_IAS = { 0.10151723027229, 0.12356453388929, 0.097174972295761,},
        L_IPS = { -0.10142460465431, 0.03229783475399, 0.11548238992691,},
        R_IAS = { 0.10771892219782, -0.12875206768513, 0.091356910765171,},
        R_IPS = { -0.1008684784174, -0.039291143417358, 0.11349472403526,},
      },
      name = "Pelvis",
      parent = "ROOT",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.245, 0.30625, 0.248713,},
          mesh_center = { 0, 0, 0.041452,},
          src = "pelvis.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.121595,},
        inertia = {
          { 0.120875, 0, 0,},
          { 0, 0.117622, 0,},
          { 0, 0, 0.023298,},
        },
        mass = 7.8334,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, -0.093, 0,},
      },
      markers = {
        R_FLE = { 0.0069601032882929, -0.042251028120518, -0.352299451828,},
        R_FME = { -0.0068211993202567, 0.075562186539173, -0.36293596029282,},
        R_FTC = { 0.012683437205851, -0.078521303832531, -0.11298349499702,},
      },
      name = "Thigh_R",
      parent = "Pelvis",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.168321, 0.134656, 0.403969,},
          mesh_center = { 0, 0, -0.168321,},
          src = "thighR.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.174376,},
        inertia = {
          { 0.029177, 0, 0,},
          { 0, 0.028309, 0,},
          { 0, 0, 0.003464,},
        },
        mass = 2.5493,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.336641,},
      },
      markers = {
        R_FAL = { 0.007205315399915, -0.04474076256156, -0.41474971175194,},
        R_FAX = { 0.0062393802218139, -0.04731772840023, -0.082377709448338,},
        R_TAM = { 0.012453241273761, 0.038268443197012, -0.40135180950165,},
        R_TTC = { 0.041200470179319, 0.012856902554631, -0.09606297314167,},
      },
      name = "Shank_R",
      parent = "Thigh_R",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.120204, 0.120204, 0.480817,},
          mesh_center = { 0, 0, -0.20034,},
          src = "shankR.obj",
        },
      },
    },
    {
      body = {
        com = { 0.002717, 0, -0.0375,},
        inertia = {
          { 0.002659, 0, 0,},
          { 0, 0.002315, 0,},
          { 0, 0, 0.000575,},
        },
        mass = 0.6837,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.400681,},
      },
      markers = {
        R_FCC = { -0.049563121050596, -0.0013734656386077, -0.058902841061354,},
        R_FM1 = { 0.12663945555687, 0.047132275998592, -0.026738381013274,},
        R_FM2 = { 0.17162348330021, -0.015348057262599, -0.014834091067314,},
        R_FM5 = { 0.10802364349365, -0.05714114382863, -0.066827692091465,},
      },
      name = "Foot_R",
      parent = "Shank_R",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.208562, 0.113, 0.075,},
          mesh_center = { 0.023281, 0, -0.0375,},
          src = "footR.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.121595,},
        inertia = {
          { 0.120875, 0, 0,},
          { 0, 0.117622, 0,},
          { 0, 0, 0.023298,},
        },
        mass = 7.8334,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0.093, 0,},
      },
      markers = {
        L_FLE = { 0.014648776501417, 0.039603006094694, -0.35856136679649,},
        L_FME = { -0.01824190095067, -0.073185376822948, -0.3641149699688,},
        L_FTC = { 0.02827356941998, 0.073947221040726, -0.11764673888683,},
      },
      name = "Thigh_L",
      parent = "Pelvis",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.168321, 0.134656, 0.403969,},
          mesh_center = { 0, 0, -0.168321,},
          src = "thighL.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.174376,},
        inertia = {
          { 0.029177, 0, 0,},
          { 0, 0.028309, 0,},
          { 0, 0, 0.003464,},
        },
        mass = 2.5493,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.336641,},
      },
      markers = {
        L_FAL = { 0.0023673840332776, 0.043155077844858, -0.41253766417503,},
        L_FAX = { 0.0041284533217549, 0.053544659167528, -0.077422261238098,},
        L_TAM = { 0.013843146152794, -0.040283098816872, -0.40212842822075,},
        L_TTC = { 0.039655886590481, -0.0074178157374263, -0.097223445773125,},
      },
      name = "Shank_L",
      parent = "Thigh_L",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.120204, 0.120204, 0.480817,},
          mesh_center = { 0, 0, -0.20034,},
          src = "shankL.obj",
        },
      },
    },
    {
      body = {
        com = { 0.002717, 0, -0.0375,},
        inertia = {
          { 0.002659, 0, 0,},
          { 0, 0.002315, 0,},
          { 0, 0, 0.000575,},
        },
        mass = 0.6837,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.400681,},
      },
      markers = {
        L_FCC = { -0.049438431859016, -0.0033535058610141, -0.061509523540735,},
        L_FM1 = { 0.13220670819283, -0.043939426541328, -0.026287345215678,},
        L_FM2 = { 0.17481558024883, 0.015778796747327, -0.018527276813984,},
        L_FM5 = { 0.10409911721945, 0.057524062693119, -0.065302357077599,},
      },
      name = "Foot_L",
      parent = "Shank_L",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.208562, 0.113, 0.075,},
          mesh_center = { 0.023281, 0, -0.0375,},
          src = "footL.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, 0.102928,},
        inertia = {
          { 0.051207, 0, 0,},
          { 0, 0.034226, 0,},
          { 0, 0, 0.047038,},
        },
        mass = 7.7645,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, 0.165808,},
      },
      markers = {
        LV1 = { -0.10315482318401, -0.002230258192867, 0.02535830065608,},
        LV3 = { -0.095416620373726, -0.0056022624485195, -0.016491221264005,},
        MAI = { -0.10759778320789, 0.0016531756846234, 0.2332865446806,},
        SXS = { 0.085490122437477, -0.006412336602807, 0.20035117864609,},
      },
      name = "MiddleTrunk",
      parent = "Pelvis",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.2205, 0.2695, 0.281326,},
          mesh_center = { 0, 0, 0.046888,},
          src = "middleTrunk.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, 0.103103,},
        inertia = {
          { 0.077145, 0, 0,},
          { 0, 0.035026, 0,},
          { 0, 0, 0.071619,},
        },
        mass = 8.1885,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, 0.187551,},
      },
      markers = {
        CV7 = { -0.064132116734982, 0.0034495270811021, 0.17492783069611,},
        SJN = { 0.073388375341892, -0.0057885749265552, 0.091342225670815,},
        TV2 = { -0.09244479984045, 0.001520021003671, 0.14115063846111,},
      },
      name = "UpperTrunk",
      parent = "MiddleTrunk",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.16695, 0.318, 0.218703,},
          mesh_center = { 0, 0, 0.083315,},
          src = "upperTrunk.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, 0.114855,},
        inertia = {
          { 0.012887, 0, 0,},
          { 0, 0.015271, 0,},
          { 0, 0, 0.011954,},
        },
        mass = 3.5404,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, 0.208288,},
      },
      markers = {
        L_HEAD = { -0.00026033207541332, 0.08052821457386, 0.14027436077595,},
        R_HEAD = { 0.0011763107031584, -0.080123014748096, 0.1429828107357,},
        SGL = { 0.092376530170441, -0.0017631795490161, 0.18659326434135,},
      },
      name = "Head",
      parent = "UpperTrunk",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.178105, 0.178105, 0.233762,},
          mesh_center = { 0, 0, 0.089052,},
          src = "head.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.144607,},
        inertia = {
          { 0.006597, 0, 0,},
          { 0, 0.00577, 0,},
          { 0, 0, 0.00187,},
        },
        mass = 1.3515,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, -0.159, 0.156288,},
      },
      markers = {
        R_HLE = { 0.041800994426012, -0.028572579845786, -0.2492528706789,},
        R_HUM = { -0.015371179208159, -0.031719449907541, -0.10547573119402,},
        R_SAE = { -0.027647502720356, -0.033248860388994, 0.028541136533022,},
      },
      name = "UpperArm_R",
      parent = "UpperTrunk",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.125658, 0.100526, 0.276448,},
          mesh_center = { 0, 0, -0.125658,},
          src = "upperArmR.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.110077,},
        inertia = {
          { 0.002905, 0, 0,},
          { 0, 0.002816, 0,},
          { 0, 0, 0.000377,},
        },
        mass = 0.7314,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.251316,},
      },
      markers = {
        R_RSP = { -0.02734780870378, 0.0014611760852858, -0.2371539324522,},
        R_USP = { 0.036448087543249, 0.00054786197142676, -0.22659465670586,},
      },
      name = "LowerArm_R",
      parent = "UpperArm_R",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.09658, 0.072435, 0.265595,},
          mesh_center = { 0, 0, -0.120725,},
          src = "lowerArmR.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.053254,},
        inertia = {
          { 0.000427, 0, 0,},
          { 0, 0.00031, 0,},
          { 0, 0, 0.00017,},
        },
        mass = 0.2968,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.24145,},
      },
      markers = {
        R_HM2 = { 3.6847512092208e-05, 6.0907848819625e-05, -0.089718461036682,},
      },
      name = "Hand_R",
      parent = "LowerArm_R",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.108776, 0.046618, 0.155394,},
          mesh_center = { 0, 0, -0.077697,},
          src = "handR.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.144607,},
        inertia = {
          { 0.006597, 0, 0,},
          { 0, 0.00577, 0,},
          { 0, 0, 0.00187,},
        },
        mass = 1.3515,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0.159, 0.156288,},
      },
      markers = {
        L_HLE = { 0.043509311974049, 0.02375622279942, -0.22889351844788,},
        L_HUM = { -0.0017371923895553, 0.028008013963699, -0.097883343696594,},
        L_SAE = { -0.012175195850432, 0.027587383985519, 0.045671228319407,},
      },
      name = "UpperArm_L",
      parent = "UpperTrunk",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.125658, 0.100526, 0.276448,},
          mesh_center = { 0, 0, -0.125658,},
          src = "upperArmL.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.110077,},
        inertia = {
          { 0.002905, 0, 0,},
          { 0, 0.002816, 0,},
          { 0, 0, 0.000377,},
        },
        mass = 0.7314,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.251316,},
      },
      markers = {
        L_RSP = { 0.021635254845023, -0.0060536777600646, -0.23349913954735,},
        L_USP = { -0.042269576340914, -0.0049579245969653, -0.2239653468132,},
      },
      name = "LowerArm_L",
      parent = "UpperArm_L",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.09658, 0.072435, 0.265595,},
          mesh_center = { 0, 0, -0.120725,},
          src = "lowerArmL.obj",
        },
      },
    },
    {
      body = {
        com = { 0, 0, -0.053254,},
        inertia = {
          { 0.000427, 0, 0,},
          { 0, 0.00031, 0,},
          { 0, 0, 0.00017,},
        },
        mass = 0.2968,
      },
      joint = {
        { 0, 1, 0, 0, 0, 0,},
        { 1, 0, 0, 0, 0, 0,},
        { 0, 0, 1, 0, 0, 0,},
      },
      joint_frame = {
        E = {
          { 1, 0, 0,},
          { 0, 1, 0,},
          { 0, 0, 1,},
        },
        r = { 0, 0, -0.24145,},
      },
      markers = {
        L_HM2 = { 0.00013292982475832, -5.9169502492296e-05, -0.08435532450676,},
      },
      name = "Hand_L",
      parent = "LowerArm_L",
      visuals = {
        {
          color = { 0.2, 0.7, 0.3,},
          dimensions = { 0.108776, 0.046618, 0.155394,},
          mesh_center = { 0, 0, -0.077697,},
          src = "handL.obj",
        },
      },
    },
  },
  gravity = { 0, 0, -9.81,},
  metadata = {
    {
      scaling_used = { "deLeva1996_segmentedTrunk",},
      subject_age = { 65,},
      subject_footWidth = { 0.113,},
      subject_gender = { "female",},
      subject_heelAnkleXOffset = { 0.081,},
      subject_heelAnkleZOffset = { 0.075,},
      subject_height = { 1.58,},
      subject_hipCenterWidth = { 0.186,},
      subject_pelvisWidth = { 0.245,},
      subject_shoulderCenterWidth = { 0.318,},
      subject_shoulderNeckZOffset = { 0.052,},
      subject_weight = { 53,},
    },
  },
  points = {
    {
      body = "Pelvis",
      name = "Pelvis_L",
      point = { 0, 0.165808, 0.116066,},
    },
    {
      body = "Pelvis",
      name = "Pelvis_R",
      point = { 0, -0.165808, 0.116066,},
    },
    {
      body = "Pelvis",
      name = "Pelvis_Back",
      point = { -0.124356, 0, 0.116066,},
    },
    {
      body = "Pelvis",
      name = "Pelvis_Front",
      point = { 0.124356, 0, 0.116066,},
    },
    {
      body = "Thigh_R",
      name = "Thigh_R",
      point = { 0.067328, 0, -0.252481,},
    },
    {
      body = "Foot_R",
      name = "Heel_Medial_R",
      point = { -0.081, 0.0565, -0.075,},
    },
    {
      body = "Foot_R",
      name = "Heel_Lateral_R",
      point = { -0.081, -0.0565, -0.075,},
    },
    {
      body = "Foot_R",
      name = "Toe_R",
      point = { 0.075422, 0, -0.075,},
    },
    {
      body = "Thigh_L",
      name = "Thigh_L",
      point = { 0.067328, 0, -0.252481,},
    },
    {
      body = "Foot_L",
      name = "Heel_Medial_L",
      point = { -0.081, -0.0565, -0.075,},
    },
    {
      body = "Foot_L",
      name = "Heel_Lateral_L",
      point = { -0.081, 0.0565, -0.075,},
    },
    {
      body = "Foot_L",
      name = "Toe_L",
      point = { 0.075422, 0, -0.075,},
    },
    {
      body = "UpperTrunk",
      name = "UpperTrunk_Front",
      point = { 0.083315, 0, 0.062486,},
    },
    {
      body = "UpperTrunk",
      name = "UpperTrunk_Back",
      point = { -0.09373, 0, 0.083315,},
    },
    {
      body = "Hand_R",
      name = "ProximalMetacarpal_Medial_R",
      point = { -0.031079, 0.023309, -0.031079,},
    },
    {
      body = "Hand_R",
      name = "ProximalMetacarpal_Lateral_R",
      point = { 0.031079, 0.023309, -0.031079,},
    },
    {
      body = "Hand_R",
      name = "DistalMetacarpal_Medial_R",
      point = { -0.031079, 0.023309, -0.093236,},
    },
    {
      body = "Hand_R",
      name = "DistalMetacarpal_Lateral_R",
      point = { 0.031079, 0.023309, -0.093236,},
    },
    {
      body = "Hand_L",
      name = "ProximalMetacarpal_Medial_L",
      point = { -0.031079, -0.023309, -0.031079,},
    },
    {
      body = "Hand_L",
      name = "ProximalMetacarpal_Lateral_L",
      point = { 0.031079, -0.023309, -0.031079,},
    },
    {
      body = "Hand_L",
      name = "DistalMetacarpal_Medial_L",
      point = { -0.031079, -0.023309, -0.093236,},
    },
    {
      body = "Hand_L",
      name = "DistalMetacarpal_Lateral_L",
      point = { 0.031079, -0.023309, -0.093236,},
    },
  },
}