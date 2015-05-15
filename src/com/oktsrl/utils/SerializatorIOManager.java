package com.oktsrl.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class SerializatorIOManager {

	public Object loadObject(final String file) throws Exception {
		ObjectInputStream in = null;

		try {
			in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(
					new File(file))));
			final Object o = in.readObject();
			in.close();

			return o;
		} catch (final Exception ex) {
			try {
				in.close();
			} catch (final Exception ex2) {
				ex2.printStackTrace();
			}

			throw ex;
		}
	}

	public void storeObject(final Object o, final String file) throws Exception {
		ObjectOutputStream out = null;

		try {
			final GZIPOutputStream gz = new GZIPOutputStream(
					new FileOutputStream(new File(file)));
			out = new ObjectOutputStream(gz);
			out.writeObject(o);
			gz.finish();
			out.close();
		} catch (final Exception ex) {
			try {
				out.close();
			} catch (final Exception ex2) {
				ex2.printStackTrace();
			}

			throw ex;
		}
	}
}