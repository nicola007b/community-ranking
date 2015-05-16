package com.oktsrl.utils;

public final class Action {
    public int user;
    public int item;
    public double rating;
    
    public int getUser() {
        return user;
    }

    public int getItem() {
        return item;
    }

    public double getRating() {
        return rating;
    }

    @Override
    public String toString() {
        return String.format("action: user=%d, item=%d, rating=%.1f", user, item , rating);
    }

    public Action() {
    }

    public Action(int user, int item, double rating) {
        this.user = user;
        this.item = item;
        this.rating = rating;
    }
}
