return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {25.0},
	 subject_height = {1.73},
	 subject_weight = {78.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.2490},
	 subject_hipCenterWidth = {0.1890},
	 subject_shoulderCenterWidth = {0.2430},
	 subject_heelAnkleXOffset = {0.0730},
	 subject_heelAnkleZOffset = {0.0960},
	 subject_shoulderNeckZOffset = {0.0360},
	 subject_footWidth = {0.1190},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Thigh_Proximal_R", body = "Thigh_R", point = {-0.073488, 0.000000, 0.036744,},},
	{name = "Thigh_Distal_R", body = "Thigh_R", point = {-0.073488, 0.000000, -0.091860,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.073000, 0.041650, -0.096000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.073000, -0.041650, -0.096000,},},
	{name = "ForeFoot_Medial_R", body = "Foot_R", point = {0.188788, 0.053550, -0.096000,},},
	{name = "ForeFoot_Lateral_R", body = "Foot_R", point = {0.188788, -0.053550, -0.096000,},},
	{name = "Thigh_Proximal_L", body = "Thigh_L", point = {-0.073488, 0.000000, 0.036744,},},
	{name = "Thigh_Distal_L", body = "Thigh_L", point = {-0.073488, 0.000000, -0.091860,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.073000, -0.041650, -0.096000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.073000, 0.041650, -0.096000,},},
	{name = "ForeFoot_Medial_L", body = "Foot_L", point = {0.188788, -0.053550, -0.096000,},},
	{name = "ForeFoot_Lateral_L", body = "Foot_L", point = {0.188788, 0.053550, -0.096000,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.033922, 0.008480, -0.084805,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.033922, 0.008480, -0.084805,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.033922, -0.008480, -0.084805,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.033922, -0.008480, -0.084805,},},
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
		mass   = 9.726600,
		com = 	{ 0.000000, 0.000000, 0.091936,},
		inertia = 
			{{ 0.059729, 0.000000, 0.000000,},
			{ 0.000000, 0.051483, 0.000000,},
			{ 0.000000, 0.000000, 0.062802,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.135733, 0.090488, 0.108586,},
		R_IAS = 	{ 0.135733, -0.090488, 0.108586,},
		L_IPS = 	{ -0.126684, 0.054293, 0.126684,},
		R_IPS = 	{ -0.126684, -0.054293, 0.126684,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.224100, 0.249000, 0.271465,},
		mesh_center = 	{ 0.000000, 0.000000, 0.045244,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.094500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 11.528400,
		com = 	{ 0.000000, 0.000000, -0.132719,},
		inertia = 
			{{ 0.211929, 0.000000, 0.000000,},
			{ 0.000000, 0.206224, 0.000000,},
			{ 0.000000, 0.000000, 0.040848,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.073488, -0.121255,},
		R_FLE = 	{ 0.000000, -0.055116, -0.367438,},
		R_FME = 	{ 0.000000, 0.055116, -0.367438,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.146975, 0.146975, 0.440926,},
		mesh_center = 	{ 0.000000, 0.000000, -0.183719,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.367438,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.751800,
		com = 	{ 0.000000, 0.000000, -0.190329,},
		inertia = 
			{{ 0.051156, 0.000000, 0.000000,},
			{ 0.000000, 0.049634, 0.000000,},
			{ 0.000000, 0.000000, 0.006074,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.087467, -0.021867,},
		R_TTC = 	{ 0.021867, 0.000000, -0.065600,},
		R_FAL = 	{ 0.000000, -0.065600, -0.437336,},
		R_TAM = 	{ 0.000000, 0.065600, -0.437336,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.131201, 0.131201, 0.524803,},
		mesh_center = 	{ 0.000000, 0.000000, -0.218668,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.437336,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.006200,
		com = 	{ 0.018376, 0.000000, -0.048000,},
		inertia = 
			{{ 0.004662, 0.000000, 0.000000,},
			{ 0.000000, 0.004059, 0.000000,},
			{ 0.000000, 0.000000, 0.001007,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.159349, 0.068293, -0.068293,},
		R_FM2 = 	{ 0.227642, 0.000000, -0.056911,},
		R_FM5 = 	{ 0.159349, -0.068293, -0.068293,},
		R_FCC = 	{ -0.068293, 0.000000, -0.022764,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.295935, 0.119000, 0.072000,},
		mesh_center = 	{ 0.074967, 0.000000, -0.060000,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.094500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 11.528400,
		com = 	{ 0.000000, 0.000000, -0.132719,},
		inertia = 
			{{ 0.211929, 0.000000, 0.000000,},
			{ 0.000000, 0.206224, 0.000000,},
			{ 0.000000, 0.000000, 0.040848,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.073488, -0.121255,},
		L_FLE = 	{ 0.000000, 0.055116, -0.367438,},
		L_FME = 	{ 0.000000, -0.055116, -0.367438,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.146975, 0.146975, 0.440926,},
		mesh_center = 	{ 0.000000, 0.000000, -0.183719,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.367438,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.751800,
		com = 	{ 0.000000, 0.000000, -0.190329,},
		inertia = 
			{{ 0.051156, 0.000000, 0.000000,},
			{ 0.000000, 0.049634, 0.000000,},
			{ 0.000000, 0.000000, 0.006074,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.087467, -0.021867,},
		L_TTC = 	{ 0.021867, 0.000000, -0.065600,},
		L_FAL = 	{ 0.000000, 0.065600, -0.437336,},
		L_TAM = 	{ 0.000000, -0.065600, -0.437336,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.131201, 0.131201, 0.524803,},
		mesh_center = 	{ 0.000000, 0.000000, -0.218668,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.437336,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.006200,
		com = 	{ 0.018376, 0.000000, -0.048000,},
		inertia = 
			{{ 0.004662, 0.000000, 0.000000,},
			{ 0.000000, 0.004059, 0.000000,},
			{ 0.000000, 0.000000, 0.001007,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.159349, -0.068293, -0.068293,},
		L_FM2 = 	{ 0.227642, 0.000000, -0.056911,},
		L_FM5 = 	{ 0.159349, 0.068293, -0.068293,},
		L_FCC = 	{ -0.068293, 0.000000, -0.022764,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.295935, 0.119000, 0.072000,},
		mesh_center = 	{ 0.074967, 0.000000, -0.060000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.180977,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 11.427000,
		com = 	{ 0.000000, 0.000000, 0.112344,},
		inertia = 
			{{ 0.089780, 0.000000, 0.000000,},
			{ 0.000000, 0.060008, 0.000000,},
			{ 0.000000, 0.000000, 0.082471,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.120500, 0.000000, 0.102354,},
		MAI = 	{ -0.120500, 0.000000, 0.102354,},
		LV1 = 	{ -0.120500, 0.000000, 0.000000,},
		LV3 = 	{ -0.120500, 0.000000, -0.051177,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.241000, 0.249000, 0.307063,},
		mesh_center = 	{ 0.000000, 0.000000, 0.051177,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.204708,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 12.051000,
		com = 	{ 0.000000, 0.000000, 0.112535,},
		inertia = 
			{{ 0.135256, 0.000000, 0.000000,},
			{ 0.000000, 0.061411, 0.000000,},
			{ 0.000000, 0.000000, 0.125568,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.096400, 0.000000, 0.181874,},
		SJN = 	{ 0.101220, 0.000000, 0.136406,},
		TV2 = 	{ -0.096400, 0.000000, 0.136406,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.241000, 0.255150, 0.238710,},
		mesh_center = 	{ 0.000000, 0.000000, 0.090937,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.227343,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 5.210400,
		com = 	{ 0.000000, 0.000000, 0.125363,},
		inertia = 
			{{ 0.022595, 0.000000, 0.000000,},
			{ 0.000000, 0.026774, 0.000000,},
			{ 0.000000, 0.000000, 0.020958,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.072899, 0.145799,},
		R_HEAD = 	{ 0.000000, -0.072899, 0.145799,},
		SGL    = 	{ 0.072899, 0.000000, 0.145799,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.194398, 0.194398, 0.255148,},
		mesh_center = 	{ 0.000000, 0.000000, 0.097199,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.121500, 0.191343,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.989000,
		com = 	{ 0.000000, 0.000000, -0.157836,},
		inertia = 
			{{ 0.011566, 0.000000, 0.000000,},
			{ 0.000000, 0.010117, 0.000000,},
			{ 0.000000, 0.000000, 0.003278,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.027431,},
		R_HUM = 	{ 0.000000, -0.041146, -0.183786,},
		R_HLE = 	{ 0.000000, -0.041146, -0.274307,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.137154, 0.109723, 0.301738,},
		mesh_center = 	{ 0.000000, 0.000000, -0.137154,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.274307,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.076400,
		com = 	{ 0.000000, 0.000000, -0.120147,},
		inertia = 
			{{ 0.005093, 0.000000, 0.000000,},
			{ 0.000000, 0.004938, 0.000000,},
			{ 0.000000, 0.000000, 0.000661,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ -0.052708, 0.000000, -0.263538,},
		R_RSP = 	{ 0.052708, 0.000000, -0.263538,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.105415, 0.079061, 0.289892,},
		mesh_center = 	{ 0.000000, 0.000000, -0.131769,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.263538,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.436800,
		com = 	{ 0.000000, 0.000000, -0.058125,},
		inertia = 
			{{ 0.000748, 0.000000, 0.000000,},
			{ 0.000000, 0.000544, 0.000000,},
			{ 0.000000, 0.000000, 0.000298,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.135688,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.118727, 0.050883, 0.169610,},
		mesh_center = 	{ 0.000000, 0.000000, -0.084805,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.121500, 0.191343,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.989000,
		com = 	{ 0.000000, 0.000000, -0.157836,},
		inertia = 
			{{ 0.011566, 0.000000, 0.000000,},
			{ 0.000000, 0.010117, 0.000000,},
			{ 0.000000, 0.000000, 0.003278,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.027431,},
		L_HUM = 	{ 0.000000, 0.041146, -0.137154,},
		L_HLE = 	{ 0.000000, 0.041146, -0.274307,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.137154, 0.109723, 0.301738,},
		mesh_center = 	{ 0.000000, -0.000000, -0.137154,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.274307,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.076400,
		com = 	{ 0.000000, 0.000000, -0.120147,},
		inertia = 
			{{ 0.005093, 0.000000, 0.000000,},
			{ 0.000000, 0.004938, 0.000000,},
			{ 0.000000, 0.000000, 0.000661,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.052708, 0.000000, -0.263538,},
		L_RSP = 	{ 0.052708, 0.000000, -0.263538,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.105415, 0.079061, 0.289892,},
		mesh_center = 	{ 0.000000, 0.000000, -0.131769,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.263538,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.436800,
		com = 	{ 0.000000, 0.000000, -0.058125,},
		inertia = 
			{{ 0.000748, 0.000000, 0.000000,},
			{ 0.000000, 0.000544, 0.000000,},
			{ 0.000000, 0.000000, 0.000298,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.135688,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.118727, 0.050883, 0.169610,},
		mesh_center = 	{ 0.000000, 0.000000, -0.084805,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}