/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.framebot.core;

/**
 *
 * @author wangqion
 */
public class Penalty {

    public static final int gap_extension_cost = 1;
    public static final int v = 10;   // gap opening penalty
    public static final int u = -gap_extension_cost;   // extension cost
    public static final int FRAME_SHIFT_PENALTY = 10;
    public static final int w = getIndelCost(1);

    public static int getIndelCost(int k) {
        return -(gap_extension_cost * k + v);
    }
}
