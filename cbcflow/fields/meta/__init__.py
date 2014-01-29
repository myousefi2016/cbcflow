# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.

# Fields that need input
meta_fields = [
    # Inspecting parts of fields
    "SubFunction",
    "PointEval",

    # Spatial analysis of other fields
    "FlowRate",
    "L2norm",
    "H1norm",
    "H1seminorm",
    "Linfnorm",
    "DiffL2norm",
    "DiffH1norm",
    "DiffH1seminorm",
    "DomainAvg",
    "BoundaryAvg",
    "Magnitude",
    "Minimum",
    "Maximum",

    # Transient analysis of other fields
    "TimeIntegral",
    "TimeDerivative",
    "SecondTimeDerivative",
    "RunningAvg",
    "RunningMin",
    "RunningMax",
    "RunningL2norm",
    ]

