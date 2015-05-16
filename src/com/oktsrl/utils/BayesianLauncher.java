package com.oktsrl.utils;

import com.oktsrl.BuildMatrixFactoryOKT;
import com.oktsrl.MatrixFactoryOKT;
import com.oktsrl.MatrixOKT;
import com.oktsrl.common.dataset.DatasetReader;
import com.oktsrl.inferencers.BayesianInferencer;
import com.oktsrl.models.GBModel;

import java.util.ArrayList;

public final class BayesianLauncher {
    private MatrixFactoryOKT factory;

    private MatrixOKT generateY(MatrixOKT A, int nnz) {
        int count= nnz<<1;
        MatrixOKT Y= factory.sparse(A.rowsCount(), A.columnsCount());
        int row, column, nUsers= A.rowsCount();
        while (count>0) {
            row= factory.nextRandomInt(nUsers);
            column= factory.nextRandomInt(nUsers);
            if( Y.get(row, column)==0 && A.get(row, column)==0 ) {
                Y.set(row, column, 1);
                count--;
            }
        }
        return Y;
    }

    public GBModel start(Settings settings, DatasetReader dataset) throws Exception {
        System.out.println("CUCU");

        // costruisco la matrice R
        ArrayList<Action> train_vec= new ArrayList<Action>();
        dataset.start();
        Review review;
        int nUsers= Integer.MIN_VALUE;
        int nItems= Integer.MIN_VALUE;
        while(dataset.hasNext()) {
            review= (Review) dataset.next(null);
            nItems= Math.max(nItems, review.item);
            nUsers= Math.max(nUsers, review.user);
            train_vec.add(new Action(review.user, review.item, review.rating));
        }

        if( nUsers<240 )
            factory= BuildMatrixFactoryOKT.getInstance(BuildMatrixFactoryOKT.BLAS);
        else
            factory= BuildMatrixFactoryOKT.getInstance(BuildMatrixFactoryOKT.UJMP);

        MatrixOKT R= factory.sparse(nUsers, nItems);
        for(Action action: train_vec)
            R.set(action.getUser()-1, action.getItem()-1, action.getRating());

        // genero la matrice A
        MatrixOKT AA= R.mul(R.transpose(true));
        MatrixOKT A= factory.sparse(AA.rowsCount());
        double threshold= settings.getReal("socialNetworkThreshold", 0d);
        int nnz= 0;
        for(int u= 0; u<AA.rowsCount(); u++) {
            for(int v= 0; v<AA.columnsCount(); v++) {
                if( AA.get(u, v)>threshold ) {
                    A.set(u, v, 1);
                    nnz++;
                }
            }
        }
        AA= null;
        //genero la matrice Y
        MatrixOKT Y= generateY(A, nnz);
        BayesianInferencer inferencer;
        System.out.println("thread"+settings.isTruth("multithread", false));
        if( settings.isTruth("multithread", false) )
            inferencer= InferencerModelFactory.getInferences("Bayesian", settings);
        else
            inferencer= InferencerModelFactory.getInferences("BayesianParallel", settings);
        inferencer.setSocialNetworkMatrix(A);
        inferencer.setSocialFoldMatrix(Y);
        inferencer.setTrainingMatrix(R);
        return (GBModel) inferencer.runInference();
    }
}