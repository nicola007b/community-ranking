package com.oktsrl.utils;


public final class Review {
    public int user;
    public int item;
    public double rating;
    public long timestamp;

    public Review() { }

    public Review(int user, int item, double rating, long timestamp) {
        this.user = user;
        this.item = item;
        this.rating = rating;
        this.timestamp = timestamp;
    }
}
