package com.oktsrl.utils;

import com.oktsrl.Inferencer;

import java.lang.reflect.Constructor;

public final class InferencerModelFactory {

    public static <T extends Inferencer> T getInferences(String name, Settings settings) throws Exception {
        String path= InferencerModelFactory.class.getPackage().getName();
        path= path.substring(0, path.lastIndexOf('.'));
        String className= String.format("%s.inferencers.impl.%sModelInferencerImpl", path, name);
        Class clazz= Class.forName(className);
        Constructor constructor= clazz.getConstructor(Settings.class);
        return (T) constructor.newInstance(settings);
    }
}
