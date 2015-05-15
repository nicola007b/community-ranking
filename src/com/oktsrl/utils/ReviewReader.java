package com.oktsrl.utils;


public interface ReviewReader {
    public boolean nextAction(Action action) throws Exception;
    public void start() throws Exception;
    public void stop() throws Exception;
}
