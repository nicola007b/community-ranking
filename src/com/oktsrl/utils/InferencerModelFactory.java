package com.oktsrl.utils;

import java.lang.reflect.Constructor;

import com.oktsrl.Inferencer;

public final class InferencerModelFactory {

	@SuppressWarnings("unchecked")
	public static <T extends Inferencer> T getInferences(String name,
			Settings settings) throws Exception {

		String path = InferencerModelFactory.class.getPackage().getName();
		path = path.substring(0, path.lastIndexOf('.'));

		final String className = String.format(
				"%s.inferencers.impl.%sModelInferencerImpl", path, name);

		final Class<T> clazz = (Class<T>) Class.forName(className);
		final Constructor<T> constructor = clazz.getConstructor(Settings.class);

		return constructor.newInstance(settings);
	}
}