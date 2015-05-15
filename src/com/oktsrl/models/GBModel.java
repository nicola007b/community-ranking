package com.oktsrl.models;

import com.oktsrl.Model;

public interface GBModel extends Model {
    public double getRatingUI(int user, int item);
    public double getRatingUX(int user);
    public double getRatingXI(int item);

    public double getRatingUU(int user_u, int user_v);
    public double getRatingU(int user, boolean vsUser);
    public double getMeanRating();
    public void computeAbsoluteTrustValue();
}
