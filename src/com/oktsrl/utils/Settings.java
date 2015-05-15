package com.oktsrl.utils;

import java.util.Properties;

public final class Settings {
	private final Properties data;

	public Settings() {
		data = new Properties();
	}

	public Settings(final Properties data) {
		this.data = data;
	}

	public Integer getInt(final String name) {
		return getInt(name, null);
	}

	public Integer getInt(final String name, final Integer defaultValue) {
		if (data != null && data.containsKey(name)) {
			final Object value = data.get(name);
			if (value instanceof Number)
				return ((Number) value).intValue();
			try {
				return Integer.valueOf(value.toString());
			} catch (final NumberFormatException e) {
				e.printStackTrace();
			}
		}
		return defaultValue;
	}

	public Double getReal(final String name) {
		return getReal(name, null);
	}

	public Double getReal(final String name, final Double defaultValue) {
		if (data != null && data.containsKey(name)) {
			final Object value = data.get(name);
			if (value instanceof Number)
				return ((Number) value).doubleValue();
			try {
				return Double.valueOf(value.toString());
			} catch (final NumberFormatException e) {
				e.printStackTrace();
			}
		}
		return defaultValue;
	}

	public String getString(final String name) {
		return getString(name, null);
	}

	public String getString(final String name, final String defaultValue) {
		if (data != null && data.containsKey(name)) {
			final Object value = data.get(name);
			if (value instanceof String)
				return (String) value;
		}
		return defaultValue;
	}

	public Boolean isTruth(final String name, final Boolean defaultValue) {
		if (data != null && data.containsKey(name)) {
			final Object value = data.get(name);
			if (value instanceof Boolean)
				return (Boolean) value;
			try {
				return Boolean.valueOf(value.toString());
			} catch (final Exception e) {
				e.printStackTrace();
			}
		}
		return defaultValue;
	}
}