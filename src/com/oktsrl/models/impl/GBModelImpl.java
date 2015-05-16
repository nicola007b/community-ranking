package com.oktsrl.models.impl;


import com.oktsrl.MatrixOKT;
import com.oktsrl.models.GBModel;

public final class GBModelImpl implements GBModel {
    private static final long serialVersionUID = 4432232006908441346l;

    private MatrixOKT[] ThetaAll;
    private MatrixOKT[] OmegaAll;
    private double meanRating;
    private double absoluteTrustValue;

    public double getMeanRating() {
        return meanRating;
    }

    @Override
    public double getRatingUI(int user, int item) {
        double sum= 0;
        int nEpoch= OmegaAll.length;
        for(int epoch=0; epoch<nEpoch; epoch++)
            try {
                sum+= ThetaAll[epoch].rows(user).dot(OmegaAll[epoch].rows(item));
            } catch (Exception e) {
                e.printStackTrace();
            }
        double rating= (sum/(double) nEpoch);
        if(  Double.isNaN(rating) || Double.isInfinite(rating) )
            return meanRating;
        return rating+meanRating;
    }


    @Override
    public double getRatingUX(int user) {
        double sum= 0;
        for(int item= 0; item<ThetaAll[0].rowsCount(); item++)
            sum+= getRatingUI(user, item);
        double rating= (sum/(double) ThetaAll[0].rowsCount());
        if(  Double.isNaN(rating) || Double.isInfinite(rating) )
            return meanRating;
        return rating;
    }

    @Override
    public double getRatingXI(int item) {
        double sum= 0;
        for(int user= 0; user<OmegaAll[0].rowsCount(); user++)
            sum+= getRatingUI(user, item);
        double rating= (sum/(double) OmegaAll[0].rowsCount());
        if(  Double.isNaN(rating) || Double.isInfinite(rating) )
            return meanRating;
        return rating;
    }

    @Override
    public double getRatingU(int user, boolean vsUser) {
        int nUser= OmegaAll[0].rowsCount();
        double sum= 0d;
        if( vsUser )
            for(int usr=0; usr<nUser; usr++)
                sum+= getRatingUU(usr, user);
        else
            for(int usr=0; usr<nUser; usr++)
                sum+= getRatingUU(user, usr);

        double rating= sum/nUser;
        if(  Double.isNaN(rating) || Double.isInfinite(rating) )
            return 0;
        return rating;
    }


    @Override
    public double getRatingUU(int user_u, int user_v) {
        double sum= 0;
        int nEpoch=OmegaAll.length;
        for(int epoch=0; epoch<nEpoch; epoch++)
            sum+= OmegaAll[epoch].rows(user_u).dot(OmegaAll[epoch].rows(user_v));
        double rating= sum/(((double) nEpoch)*absoluteTrustValue);
        if( Double.isNaN(rating) || Double.isInfinite(rating) ) {
            return 0;
        }
        return rating;
    }

    public void computeAbsoluteTrustValue() {
        System.out.println("Calcolo dell'absoluteTrustValue in corso...");
        int nUsers= OmegaAll[0].rowsCount();
        absoluteTrustValue= 0;
        double rating;
        for(int user_u= 0; user_u<nUsers; user_u++)
            for(int user_v= 0; user_v<nUsers; user_v++) {
                rating= Math.abs(getRatingUU(user_u, user_v));
                if( !(Double.isNaN(rating) || Double.isInfinite(rating)) )
                    absoluteTrustValue= Math.max(absoluteTrustValue, rating);
            }
        if( absoluteTrustValue==0 )
            absoluteTrustValue=1;

        System.out.println("absoluteTrustValue: "+absoluteTrustValue);
    }

    public GBModelImpl(MatrixOKT[] Omega_All, MatrixOKT[] Theta_All, double meanRating) {
        OmegaAll = Omega_All;
        ThetaAll = Theta_All;
        this.meanRating= meanRating;
    }
}
