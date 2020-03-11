return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {65.0},
	 subject_height = {1.48},
	 subject_weight = {66.00},
	 subject_gender = {"male"},
	 subject_pelvisWidth = {0.2170},
	 subject_hipCenterWidth = {0.1650},
	 subject_shoulderCenterWidth = {0.2050},
	 subject_heelAnkleXOffset = {0.0990},
	 subject_heelAnkleZOffset = {0.0920},
	 subject_shoulderNeckZOffset = {0.0330},
	 subject_footWidth = {0.1110},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Pelvis_L", body = "Pelvis", point = {0.000000, 0.123439, 0.086407,},},
	{name = "Pelvis_R", body = "Pelvis", point = {0.000000, -0.123439, 0.086407,},},
	{name = "Pelvis_Back", body = "Pelvis", point = {-0.092579, 0.000000, 0.086407,},},
	{name = "Pelvis_Front", body = "Pelvis", point = {0.092579, 0.000000, 0.086407,},},
	{name = "Thigh_R", body = "Thigh_R", point = {0.071539, -0.000000, -0.268270,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.099000, 0.055500, -0.092000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.099000, -0.055500, -0.092000,},},
	{name = "Toe_R", body = "Foot_R", point = {0.064999, 0.000000, -0.092000,},},
	{name = "Thigh_L", body = "Thigh_L", point = {0.071539, 0.000000, -0.268270,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.099000, -0.055500, -0.092000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.099000, 0.055500, -0.092000,},},
	{name = "Toe_L", body = "Foot_L", point = {0.064999, 0.000000, -0.092000,},},
	{name = "UpperTrunk_Front", body = "UpperTrunk", point = {0.082044, 0.000000, 0.061533,},},
	{name = "UpperTrunk_Back", body = "UpperTrunk", point = {-0.092300, 0.000000, 0.082044,},},
	{name = "ProximalMetacarpal_Medial_R", body = "Hand_R", point = {-0.031838, 0.023879, -0.031838,},},
	{name = "ProximalMetacarpal_Lateral_R", body = "Hand_R", point = {0.031838, 0.023879, -0.031838,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.031838, 0.023879, -0.095515,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.031838, 0.023879, -0.095515,},},
	{name = "ProximalMetacarpal_Medial_L", body = "Hand_L", point = {-0.031838, -0.023879, -0.031838,},},
	{name = "ProximalMetacarpal_Lateral_L", body = "Hand_L", point = {0.031838, -0.023879, -0.031838,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.031838, -0.023879, -0.095515,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.031838, -0.023879, -0.095515,},},
},
constraint_sets = {
},
frames = {
	{name   = "Pelvis",
	parent = "ROOT",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 7.372200,
		com = 	{ 0.000000, 0.000000, 0.047956,},
		inertia = 
			{{ 0.042487, 0.000000, 0.000000,},
			{ 0.000000, 0.034104, 0.000000,},
			{ 0.000000, 0.000000, 0.038706,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.092579, 0.061720, 0.074063,},
		R_IAS = 	{ 0.092579, -0.061720, 0.074063,},
		L_IPS = 	{ -0.086407, 0.086407, 0.086407,},
		R_IPS = 	{ -0.086407, -0.086407, 0.086407,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.217000, 0.271250, 0.185159,},
		mesh_center = 	{ 0.000000, 0.000000, 0.030860,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.082500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.345600,
		com = 	{ 0.000000, 0.000000, -0.146476,},
		inertia = 
			{{ 0.129426, 0.000000, 0.000000,},
			{ 0.000000, 0.129426, 0.000000,},
			{ 0.000000, 0.000000, 0.026546,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.071539, -0.118039,},
		R_FLE = 	{ 0.000000, -0.053654, -0.357694,},
		R_FME = 	{ 0.000000, 0.053654, -0.357694,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.178847, 0.143078, 0.429233,},
		mesh_center = 	{ 0.000000, 0.000000, -0.178847,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.357694,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.857800,
		com = 	{ 0.000000, 0.000000, -0.163946,},
		inertia = 
			{{ 0.025053, 0.000000, 0.000000,},
			{ 0.000000, 0.024065, 0.000000,},
			{ 0.000000, 0.000000, 0.004137,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.074606, -0.018651,},
		R_TTC = 	{ 0.018651, 0.000000, -0.055954,},
		R_FAL = 	{ 0.000000, -0.055954, -0.373028,},
		R_TAM = 	{ 0.000000, 0.055954, -0.373028,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.111909, 0.111909, 0.447634,},
		mesh_center = 	{ 0.000000, 0.000000, -0.186514,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.373028,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.904200,
		com = 	{ -0.002459, 0.000000, -0.046000,},
		inertia = 
			{{ 0.002856, 0.000000, 0.000000,},
			{ 0.000000, 0.002595, 0.000000,},
			{ 0.000000, 0.000000, 0.000665,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.153066, 0.021867, -0.032800,},
		R_FM2 = 	{ 0.153066, 0.000000, -0.032800,},
		R_FM5 = 	{ 0.153066, -0.021867, -0.032800,},
		R_FCC = 	{ -0.065600, 0.000000, -0.032800,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.218666, 0.111000, 0.092000,},
		mesh_center = 	{ 0.010333, 0.000000, -0.046000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.082500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.345600,
		com = 	{ 0.000000, 0.000000, -0.146476,},
		inertia = 
			{{ 0.129426, 0.000000, 0.000000,},
			{ 0.000000, 0.129426, 0.000000,},
			{ 0.000000, 0.000000, 0.026546,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.071539, -0.118039,},
		L_FLE = 	{ 0.000000, 0.053654, -0.357694,},
		L_FME = 	{ 0.000000, -0.053654, -0.357694,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.178847, 0.143078, 0.429233,},
		mesh_center = 	{ 0.000000, 0.000000, -0.178847,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.357694,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.857800,
		com = 	{ 0.000000, 0.000000, -0.163946,},
		inertia = 
			{{ 0.025053, 0.000000, 0.000000,},
			{ 0.000000, 0.024065, 0.000000,},
			{ 0.000000, 0.000000, 0.004137,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.074606, -0.018651,},
		L_TTC = 	{ 0.018651, 0.000000, -0.055954,},
		L_FAL = 	{ 0.000000, 0.055954, -0.373028,},
		L_TAM = 	{ 0.000000, -0.055954, -0.373028,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.111909, 0.111909, 0.447634,},
		mesh_center = 	{ 0.000000, 0.000000, -0.186514,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.373028,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.904200,
		com = 	{ -0.002459, 0.000000, -0.046000,},
		inertia = 
			{{ 0.002856, 0.000000, 0.000000,},
			{ 0.000000, 0.002595, 0.000000,},
			{ 0.000000, 0.000000, 0.000665,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.153066, -0.021867, -0.032800,},
		L_FM2 = 	{ 0.153066, 0.000000, -0.032800,},
		L_FM5 = 	{ 0.153066, 0.021867, -0.032800,},
		L_FCC = 	{ -0.065600, 0.000000, -0.032800,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.218666, 0.111000, 0.092000,},
		mesh_center = 	{ 0.010333, 0.000000, -0.046000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.123439,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 10.777800,
		com = 	{ 0.000000, 0.000000, 0.100380,},
		inertia = 
			{{ 0.083465, 0.000000, 0.000000,},
			{ 0.000000, 0.052700, 0.000000,},
			{ 0.000000, 0.000000, 0.078687,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.091287, 0.000000, 0.091287,},
		MAI = 	{ -0.091287, 0.000000, 0.091287,},
		LV1 = 	{ -0.091287, 0.000000, 0.000000,},
		LV3 = 	{ -0.091287, 0.000000, -0.045644,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.195300, 0.238700, 0.273862,},
		mesh_center = 	{ 0.000000, 0.000000, 0.045644,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.182575,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 10.533600,
		com = 	{ 0.000000, 0.000000, 0.101202,},
		inertia = 
			{{ 0.113015, 0.000000, 0.000000,},
			{ 0.000000, 0.045379, 0.000000,},
			{ 0.000000, 0.000000, 0.095821,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.082044, 0.000000, 0.164088,},
		SJN = 	{ 0.086146, 0.000000, 0.123066,},
		TV2 = 	{ -0.082044, 0.000000, 0.123066,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.134531, 0.256250, 0.215366,},
		mesh_center = 	{ 0.000000, 0.000000, 0.082044,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.205111,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.580400,
		com = 	{ 0.000000, 0.000000, 0.102853,},
		inertia = 
			{{ 0.017809, 0.000000, 0.000000,},
			{ 0.000000, 0.019247, 0.000000,},
			{ 0.000000, 0.000000, 0.013214,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.061737, 0.123473,},
		R_HEAD = 	{ 0.000000, -0.061737, 0.123473,},
		SGL    = 	{ 0.061737, 0.000000, 0.123473,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.164631, 0.164631, 0.216078,},
		mesh_center = 	{ 0.000000, 0.000000, 0.082315,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.102500, 0.172111,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.788600,
		com = 	{ 0.000000, 0.000000, -0.137755,},
		inertia = 
			{{ 0.008275, 0.000000, 0.000000,},
			{ 0.000000, 0.007372, 0.000000,},
			{ 0.000000, 0.000000, 0.002543,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.023866,},
		R_HUM = 	{ 0.000000, -0.035799, -0.159902,},
		R_HLE = 	{ 0.000000, -0.035799, -0.238660,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.119330, 0.095464, 0.262526,},
		mesh_center = 	{ 0.000000, 0.000000, -0.119330,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.238660,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.069200,
		com = 	{ 0.000000, 0.000000, -0.104203,},
		inertia = 
			{{ 0.004227, 0.000000, 0.000000,},
			{ 0.000000, 0.003897, 0.000000,},
			{ 0.000000, 0.000000, 0.000812,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ 0.045563, 0.000000, -0.227816,},
		R_RSP = 	{ -0.045563, 0.000000, -0.227816,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.091126, 0.068345, 0.250598,},
		mesh_center = 	{ 0.000000, 0.000000, -0.113908,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.227816,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.402600,
		com = 	{ 0.000000, 0.000000, -0.057691,},
		inertia = 
			{{ 0.000846, 0.000000, 0.000000,},
			{ 0.000000, 0.000563, 0.000000,},
			{ 0.000000, 0.000000, 0.000345,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.127353,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.111434, 0.047757, 0.159192,},
		mesh_center = 	{ 0.000000, 0.000000, -0.079596,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.102500, 0.172111,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.788600,
		com = 	{ 0.000000, 0.000000, -0.137755,},
		inertia = 
			{{ 0.008275, 0.000000, 0.000000,},
			{ 0.000000, 0.007372, 0.000000,},
			{ 0.000000, 0.000000, 0.002543,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.023866,},
		L_HUM = 	{ 0.000000, 0.035799, -0.119330,},
		L_HLE = 	{ 0.000000, 0.035799, -0.238660,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.119330, 0.095464, 0.262526,},
		mesh_center = 	{ 0.000000, -0.000000, -0.119330,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.238660,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.069200,
		com = 	{ 0.000000, 0.000000, -0.104203,},
		inertia = 
			{{ 0.004227, 0.000000, 0.000000,},
			{ 0.000000, 0.003897, 0.000000,},
			{ 0.000000, 0.000000, 0.000812,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.045563, 0.000000, -0.227816,},
		L_RSP = 	{ 0.045563, 0.000000, -0.227816,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.091126, 0.068345, 0.250598,},
		mesh_center = 	{ 0.000000, 0.000000, -0.113908,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.227816,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.402600,
		com = 	{ 0.000000, 0.000000, -0.057691,},
		inertia = 
			{{ 0.000846, 0.000000, 0.000000,},
			{ 0.000000, 0.000563, 0.000000,},
			{ 0.000000, 0.000000, 0.000345,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.127353,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.111434, 0.047757, 0.159192,},
		mesh_center = 	{ 0.000000, 0.000000, -0.079596,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}